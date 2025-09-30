# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


try:
    from dask.distributed import Client, LocalCluster, Security
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.security' requires the "
        + "third-party package 'dask' as a dependency!\n"
        + "Please install this package into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/dask/\n"
    )
    raise

import collections
import logging
import os
import shutil
import subprocess

from pathlib import Path
from typing import (
    Dict,
    Generic,
    Iterable,
    List,
    NoReturn,
    Optional,
    TypeVar,
    Union,
)

from pyrosetta.distributed.cluster.hkdf import MaskedBytes, derive_task_key
from pyrosetta.distributed.cluster.serialization import NonceCache
from pyrosetta.distributed.cluster.worker_plugins import TaskSecurityPlugin


G = TypeVar("G")


class SecurityIO(Generic[G]):
    def _setup_task_security_plugin(self, clients: Dict[int, Client]) -> None:
        """Setup tast security worker plugin(s)."""
        prk = MaskedBytes(derive_task_key(os.urandom(32), self.instance_id))
        self._register_task_security_plugin(clients, prk)
        self.serializer.prk = prk
        assert self.serializer.__getstate__()["prk"] is None, "Pseudo-random key is not hidden on host serializer."
        if self.with_nonce:
            self.nonce_cache.prk = self.serializer.prk
            # The `NonceCache.__getstate__` method must be overridden after running `self._register_task_security_plugin()`
            # since WorkerPlugin registration first pickles `NonceCache.prk` with the default `__getstate__` method:
            self.nonce_cache.__getstate__ = NonceCache._get_state.__get__(self.nonce_cache, self.nonce_cache.__class__)
            assert self.nonce_cache.__getstate__()["prk"] is None, "Pseudo-random key is not hidden on host nonce cache."

    def _register_task_security_plugin(self, clients: Dict[int, Client], prk: MaskedBytes) -> None:
        """Register `TaskSecurityPlugin` as a dask worker plugin on dask clients."""
        for client in clients.values():
            plugin = TaskSecurityPlugin(self.instance_id, prk, self.max_nonce)
            plugin.idempotent = True # Never re-register plugin
            if hasattr(client, "register_plugin"):
                client.register_plugin(plugin=plugin, name=self.instance_id)
            else: # Deprecated since dask version 2023.9.2
                client.register_worker_plugin(plugin=plugin, name=self.instance_id, nanny=False)

    def _clients_dict_has_security(self) -> bool:
        """
        Test if the `self.clients_dict` has security enabled on all clients, excluding
        clients with `LocalCluster` clusters.
        """
        assert len(self.clients_dict) > 0, "No clients in `self.clients_dict` to test for `security` attribute."
        for _client in self.clients_dict.values():
            if not isinstance(_client.cluster, LocalCluster) and not isinstance(_client.security, Security):
                _has_security = False
                break
        else:
            _has_security = True

        return _has_security

    def _setup_with_nonce(self) -> bool:
        """Post-init hook to setup the `PyRosettaCluster().with_nonce` instance attribute."""
        if self.clients_dict:
            with_nonce = not self._clients_dict_has_security()
            if with_nonce:
                logging.warning(
                    "A PyRosettaCluster input client with a remote cluster does not have security enabled, "
                    + "so PyRosettaCluster is enabling nonce caching for security on the host and all worker "
                    + "processes. Please set the `PyRosettaCluster(max_nonce=...)` keyword argument parameter "
                    + "to limit nonce cache memory usage. To silence this warning, enable dask security on all "
                    + "PyRosettaCluster input clients that are not instances of `dask.distributed.LocalCluster`."
                )
        else:
            if bool(self.scheduler):
                with_nonce = not bool(self.security)
            else:
                with_nonce = False

        return with_nonce


def generate_dask_tls_security(
    output_dir: str = ".",
    common_name: str = "dask_tls_security",
    days: int = 365,
    openssl_bin: str = "openssl",
    overwrite: bool = False,
    san_dns: Optional[Iterable[str]] = None,
    san_ip: Optional[Iterable[str]] = None,
    cleanup: bool = True,
) -> Union[Security, NoReturn]:
    """
    Create cryptographic certificates and private keys for securing a dask cluster,
    and return a dask `Security` object that can be passed directly to the
    `PyRosettaCluster(security=...)` keyword argument parameter.

    This function uses the `openssl` command-line tool to generate the following:

    - A "certificate authority" certificate and key:
      These act as a trusted "parent" identity used to sign other certificates.
      Files: `ca.pem` (certificate), `ca.key` (private key).

    - A "leaf" certificate and key:
      These represent the actual dask processes (scheduler, workers, and client).
      Files: `tls.crt` (certificate), `tls.key` (private key).

    By default, the leaf certificate is signed by the certificate authority,
    meaning that any process configured with this authority will trust the leaf
    certificate as valid.

    All generated files are placed in the `output_dir` keyword argument parameter,
    which defaults to the current working directory.

    Keyword Args:
        output_dir: A `str` object representing the directory where all certificate and
            key files will be written. The directory will be created if it does not exist.
            All generated files (CA certificate, leaf certificate, leaf private key, and
            optional bookkeeping files) are output to this single directory. Therefore,
            for a distributed dask setup, this directory must be readable by the scheduler,
            workers, and client processes, either via a shared filesystem or via copying and
            mounting (e.g., if using Docker, Apptainer, or other container applications).
            Default: "."
        common_name: A `str` object representing the "Common Name" placed inside the leaf
            certificate. This is a human-readable identifier that typically names the
            system or service to which the certificate belongs.
            Default: "dask_tls_security"
        days: An `int` object representing the number of days the certificates will be
            valid before expiring.
            Default: 365
        openssl_bin: A `str` object representing the path or name of the `openssl`
            executable. If the OpenSSL executable is not in the system "PATH"
            environment variable, then the full path must be provided.
            Default: "openssl"
        overwrite: A `bool` object specifying whether or not to overwrite existing
            files in 'output_dir' keyword argument parameter. If `True` is provided,
            the same filenames will be deleted and replaced with newly generated ones.
            If `False` is provided, then existing files are re-used.
            Default: False
        san_dns: An optional iterable of `str` object representing a list of hostnames
            (e.g., `["localhost", "cluster.example.com"]`) that should be accepted when
            verifying the certificate. These are included in an extension field called
            "Subject Alternative Names".
            Default: None
        san_ip: An optional iterable of `str` object representing a list of IP addresses
            (e.g., `["127.0.0.1", "111.111.111.1"]`) that should be accepted when verifying
            the certificate. These are also included in the "Subject Alternative Names" field.
            Default: None
        cleanup: An optional `bool` object specifying whether or not to delete the `index.txt`
            and `serial` bookkeeping files used by OpenSSL.
            Default: True

    Returns:
        A `dask.distributed.Security()` instance configured to require encryption (with
        `require_encryption=True`) and configured to use the generated certificates
        and keys for the scheduler, workers, and client.


    Examples:
        Generate a new set of certificates and a configured dask `Security` object:
        ```
        security = pyrosetta.distributed.cluster.generate_dask_tls_security(
            output_dir="./dask_certs",
            common_name="my-cluster",
            san_dns=["localhost", "my-host.local"],
            san_ip=["127.0.0.1"],
            cleanup=False,
        )
        ```
        After running this function, the directory `./dask_certs` will contain:
        - `ca.pem`: certificate authority certificate (used by dask)
        - `ca.key`: certificate authority private key
        - `tls.crt`: leaf certificate (used by dask)
        - `tls.key`: leaf private key (used by dask)
        - `index.txt`, `serial`, and `ca.cnf`: bookkeeping files used by OpenSSL (with `cleanup=False`)

        Then use the configured dask `Security` object with PyRosettaCluster:
        ```
        PyRosettaCluster(security=security, ...).distribute(...)
        ```

    Additional Notes:
        - A "certificate authority" (CA) act as a trusted "parent" identity that
        confirms whether a certificate is real. In this function, the user generates
        their own local CA for the dask cluster.
        - A "leaf certificate" is the actual identity used by a running process
        (i.e., the scheduler, a worker, or a client).
        - "Subject Alternative Names" (SANs) are extra hostnames or IP addresses
        for which the certificate is valid. This enables the user to connect using
        either a machine name or an IP address without validation errors.
        - File permissions are automatically set for private keys using `chmod 600`
        so they are restricted to the owner (read/write only) for basic security.
        - This function generates all necessary files in a single directory. For
        proper TLS validation in a distributed dask setup, the CA certificate must
        be accessible from all nodes (i.e., the scheduler, workers, and client).
        Leaf certificates and keys must be accessible by the process using them.
        For example, all files can be placed in a common directory from which all
        processes can read, or the directory can be mounted (e.g., if using Docker,
        Apptainer, or other container applications).
        - If `cleanup=False` and the same directory is used for multiple function
        calls, then OpenSSL may create additional files in the output directory
        (e.g., `*.pem`, `index.txt.attr`, `index.txt.old`, and `serial.old`).
        These are simply bookkeeping files used internally by OpenSSL and are not
        required by dask, so they can be safely deleted after the leaf certificate
        has been issued.
    """
    def _run(cmd: List[str]) -> Optional[NoReturn]:
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as ex:
            raise RuntimeError(
                "OpenSSL command failed!\n"
                + f"Command: {' '.join(cmd) if isinstance(cmd, list) else cmd}\n"
                + f"Return code: {ex.returncode}\n"
                + f"STDOUT:\n{ex.stdout}\n"
                + f"STDERR:\n{ex.stderr}\n"
            ) from ex

    def _chmod_600(path: Path) -> None:
        try:
            os.chmod(path, 0o600)
        except Exception:
            pass # Ignored on platforms that do not support `chmod`

    def _maybe_unlink(path: Path) -> None:
        try:
            if path.exists():
                path.unlink()
        except Exception:
            pass # Ignore if files cannot be deleted

    # Check types
    if not isinstance(output_dir, str):
        raise TypeError(
            "The 'output_dir' keyword argument parameter must be of type `str`. "
            + f"Received: {type(output_dir)}"
        )
    if not isinstance(common_name, str):
        raise TypeError(
            "The 'common_name' keyword argument parameter must be of type `str`. "
            f"Received: {type(common_name)}"
        )
    if not isinstance(days, int):
        raise TypeError(
            "The 'days' keyword argument parameter must be of type `int`. "
            + f"Received: {type(days)}"
        )
    if not isinstance(openssl_bin, str):
        raise TypeError(
            "The 'openssl_bin' keyword argument parameter must be of type `str`. "
            + f"Received: {type(openssl_bin)}"
        )
    if not isinstance(overwrite, bool):
        raise TypeError(
            "The 'overwrite' keyword argument parameter must be of type `bool`. "
            + f"Received: {type(overwrite)}"
        )
    if san_dns is not None:
        if not isinstance(san_dns, collections.abc.Iterable):
            raise TypeError(
                "The 'san_dns' keyword argument parameter must be an iterable. "
                + f"Received: {type(san_dns)}"
            )
        else:
            for obj in san_dns:
                if not isinstance(obj, str):
                    raise TypeError(
                        "The 'san_dns' keyword argument parameter must be an iterable of `str` objects. "
                        + f"Received: {type(obj)}"
                    )
    if san_ip is not None:
        if not isinstance(san_ip, collections.abc.Iterable):
            raise TypeError(
                "The 'san_ip' keyword argument parameter must be an iterable. "
                + f"Received: {type(san_ip)}"
            )
        else:
            for obj in san_ip:
                if not isinstance(obj, str):
                    raise TypeError(
                        "The 'san_ip' keyword argument parameter must be an iterable of `str` objects. "
                        + f"Received: {type(obj)}"
                    )
    if not isinstance(cleanup, bool):
        raise TypeError(
            "The 'cleanup' keyword argument parameter must be of type `bool`. "
            + f"Received: {type(cleanup)}"
        )

    openssl = shutil.which(openssl_bin)
    if not openssl:
        raise FileNotFoundError(
            f"Could not find the provided OpenSSL executable. Please install OpenSSL or pass the full "
            + f"path to the 'openssl_bin' keyword argument parameter. Received: '{openssl_bin}'"
        )

    outdir = Path(os.path.expanduser(output_dir)).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    ca_cert = outdir / "ca.pem"
    ca_key = outdir / "ca.key"
    leaf_crt = outdir / "tls.crt"
    leaf_key = outdir / "tls.key"
    leaf_csr = outdir / "tls.csr"
    leaf_ext = outdir / "leaf_ext.cnf"

    if overwrite:
        for path in (ca_cert, ca_key, leaf_crt, leaf_key, leaf_csr, leaf_ext):
            _maybe_unlink(path)

    if ca_key.exists():
        _chmod_600(ca_key)

    # Create CA if it does not exist
    if not (ca_cert.exists() and ca_key.exists()):
        _run(
            [
                openssl,
                "req",
                "-x509",
                "-newkey", "rsa:2048",
                "-days", str(days),
                "-nodes",
                "-keyout", str(ca_key),
                "-out", str(ca_cert),
                "-verbose",
                "-subj", "/CN=dask_security_ca",
            ]
        )
        _chmod_600(ca_key)

    # Create leaf key + CSR if they do not exist
    if not (leaf_crt.exists() and leaf_key.exists()):
        _run(
            [
                openssl,
                "req",
                "-new",
                "-newkey", "rsa:2048",
                "-nodes",
                "-keyout", str(leaf_key),
                "-out", str(leaf_csr),
                "-subj", f"/CN={common_name}",
                "-verbose",
            ]
        )
        _chmod_600(leaf_key)

        # Optional SANs: build extfile if requested
        san_entries = []
        if san_dns:
            san_entries += [f"DNS:{d}" for d in san_dns]
        if san_ip:
            san_entries += [f"IP:{ip}" for ip in san_ip]
        if san_entries:
            text = "[v3_req]\nsubjectAltName = " + ", ".join(san_entries) + "\n"
            leaf_ext.write_text(text)
        else:
            # Ensure extfile does not exist if no SANs requested
            if leaf_ext.exists():
                _maybe_unlink(leaf_ext)

        # Setup CA database
        index = outdir / "index.txt"
        serial = outdir / "serial"
        if not index.exists():
            index.write_text("")
        if not serial.exists():
            serial.write_text("1000\n")

        # Build a custom CA config
        ca_config = outdir / "ca.cnf"
        config_text = f"""
        [ ca ]
        default_ca = CA_default

        [ CA_default ]
        dir               = {outdir}
        new_certs_dir     = {outdir}
        database          = {index}
        serial            = {serial}
        certificate       = {ca_cert}
        private_key       = {ca_key}
        default_md        = sha256
        policy            = policy_any
        x509_extensions   = v3_req

        [ policy_any ]
        commonName        = supplied

        [ v3_req ]
        """
        if san_entries:
            config_text += "subjectAltName = " + ", ".join(san_entries) + "\n"

        ca_config.write_text(config_text)

        # Run `openssl ca` with custom configuration
        cmd = [
            openssl,
            "ca",
            "-config", str(ca_config),
            "-in", str(leaf_csr),
            "-out", str(leaf_crt),
            "-days", str(days),
            "-batch",
        ]
        _run(cmd)

        # Update permissions on produced files
        _chmod_600(leaf_key)
        if leaf_crt.exists():
            _chmod_600(leaf_crt)
        # Clean up intermediate files
        _maybe_unlink(leaf_csr)
        _maybe_unlink(leaf_ext)
        if cleanup:
            # Clean up OpenSSL bookkeeping files
            for path in (index, serial, ca_config):
                _maybe_unlink(path)
            for path in outdir.glob("*.old"):
                _maybe_unlink(path)
            for path in outdir.glob("*.attr"):
                _maybe_unlink(path)
            for path in outdir.glob("[0-9]*.pem"):
                _maybe_unlink(path)

    # Return dask `Security` object configured to use the CA and the leaf certs
    return Security(
        tls_ca_file=str(ca_cert),
        tls_client_cert=str(leaf_crt),
        tls_client_key=str(leaf_key),
        tls_scheduler_cert=str(leaf_crt),
        tls_scheduler_key=str(leaf_key),
        tls_worker_cert=str(leaf_crt),
        tls_worker_key=str(leaf_key),
        require_encryption=True,
    )
