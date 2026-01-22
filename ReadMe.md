# Notes to myself for the plastid genome structure project 🧬🌷🌿

## Necessary modules for the pipeline to run:
* Bio
* os
* time
* certifi
* argparse

## Potential errors

fetchGenbank,py
`Failed for NC_068711: <urlopen error [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: self signed certificate in certificate chain (_ssl.c:992)>`
Had to put "import certifi" and "os.environ["SSL_CERT_FILE"] = certifi.where()" for this to get fixed.