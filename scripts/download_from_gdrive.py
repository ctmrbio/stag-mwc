#!/usr/bin/env python
# Download files from shared Google Drive links
# https://stackoverflow.com/questions/25010369/wget-curl-large-file-from-google-drive
# UI improvements by Fredrik Boulund 2017

from sys import argv, exit
import argparse
import requests


def parse_args():
    desc = "Download files from Google Drive shared links."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("drive_file_id", metavar="GDRIVE_ID",
            help="File ID from Google Drive shared link.")
    parser.add_argument("-o", "--outfile", metavar="FILENAME",
            default="GDrive_download",
            help="Save download to FILENAME [%(default)s].")
    if len(argv) < 2:
        parser.print_help()
        exit(1)

    return parser.parse_args()


def download_file_from_google_drive(gdrive_id, destination):
    def get_confirm_token(response):
        for key, value in response.cookies.items():
            if key.startswith("download_warning"):
                return value
        return None

    def save_response_content(response, destination):
        CHUNK_SIZE = 32768
        with open(destination, "wb") as f:
            print("Downloading to:", destination)
            for chunk in response.iter_content(CHUNK_SIZE):
                if chunk: # filter out keep-alive new chunks
                    f.write(chunk)

    base_URL = "https://docs.google.com/uc?export=download"

    session = requests.Session()
    response = session.get(base_URL, params={"id": gdrive_id}, stream=True)
    token = get_confirm_token(response)

    if token:
        params = {"id": gdrive_id, "confirm": token}
        response = session.get(base_URL, params=params, stream=True)

    save_response_content(response, destination)


if __name__ == "__main__":
    options = parse_args()
    download_file_from_google_drive(options.drive_file_id, options.outfile)
