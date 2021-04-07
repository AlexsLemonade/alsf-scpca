#!/usr/bin/env python3
"""
Script to check md5 of files on S3.

"""

import argparse
import boto3
import hashlib

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--md5-file',
                        dest = 'md5_file')
    parser.add_argument('-p', '--s3_prefix',
                        dest = 'prefix',
                        help = 'The s3 location of the files, including bucket')
    args = parser.parse_args()

    ## Read md5 file
    md5_pairs = []
    with open(args.md5_file) as f:
        md5_pairs = [line.strip().split() for line in f]
    
    md5_hash = hashlib.md5()
    
    # set up s3
    s3 = boto3.client('s3')
    bucket, prefix = args.prefix.split('/', 1)

    for md5, filename in md5_pairs:
        md5_hash = hashlib.md5()
        # get the object and stream in chunks
        try:
            obj = s3.get_object(Bucket = bucket, Key = f"{prefix}/{filename}")
        except:
            print(f"{prefix}/{filename}: Missing")
            continue

        chunks = obj['Body'].iter_chunks(65536)
        for chunk in chunks:
            md5_hash.update(chunk)
        if md5 == md5_hash.hexdigest():
            print(f"{prefix}/{filename}: OK")
        else: 
            print(f"{prefix}/{filename}: MD5 mismatch")

if __name__ == "__main__":
    main()
