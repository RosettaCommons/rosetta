#!/usr/bin/python3

import os, os.path, sys, zlib, json, copy, subprocess, shutil, urllib.parse

in_html = sys.argv[1]
out_pdf = sys.argv[2]

wkhtmltopdf_options = {
    'enable-local-file-access': None
}

import pdfkit
#import weasyprint
#from weasyprint import HTML

pdfkit.from_file( in_html, out_pdf, options = wkhtmltopdf_options )
