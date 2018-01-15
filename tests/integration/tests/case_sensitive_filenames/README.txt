This integration test checks that you haven't inadvertantly commited files which will cause issues on case-insensitive filesystems.

It doesn't look at the current state of the directory, but just at the contents of the current commit.
