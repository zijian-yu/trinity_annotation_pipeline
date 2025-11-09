#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Download UniProt protein database
"""
import os
import gzip
import urllib.request
from tqdm import tqdm


# 自定义进度条类，用于显示下载进度
class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


# 下载文件并显示进度条
def download_file(url, filename):
    print(f"Downloading {os.path.basename(url)}...")
    with DownloadProgressBar(unit='B', unit_scale=True, miniters=1, 
                           desc=os.path.basename(url)) as t:
        urllib.request.urlretrieve(url, filename=filename, reporthook=t.update_to)
    print("Download completed.")


# 解压.gz文件
def extract_gz(gz_filename, output_filename):
    print(f"Extracting {gz_filename}...")
    with gzip.open(gz_filename, 'rb') as gz_file:
        gz_file.seek(0, 2)
        file_size = gz_file.tell()
        gz_file.seek(0)
        
        with open(output_filename, 'wb') as out_file:
            with tqdm(total=file_size, unit='B', unit_scale=True, 
                     desc=f"Extracting {output_filename}") as pbar:
                while True:
                    chunk = gz_file.read(8192)  # 8KB chunks
                    if not chunk:
                        break
                    out_file.write(chunk)
                    pbar.update(len(chunk))
    print(f"File saved as {output_filename}")


# 下载并解压UniProt Swiss-Prot数据库
def download_uniprot_pep(force_download=False, output_file="uniprots.pep"):
    """
    Download and extract UniProt Swiss-Prot database
    
    Args:
        force_download: If True, re-download even if file exists
        output_file: Output filename for the protein database
        
    Returns:
        bool: True if successful, False otherwise
    """
    url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
    gz_filename = "uniprot_sprot.fasta.gz"

    if os.path.exists(output_file) and not force_download:
        print(f"{output_file} already exists. Use force_download=True to re-download.")
        return True
    try:
        download_file(url, gz_filename)
        extract_gz(gz_filename, output_file)
        os.remove(gz_filename)
        print(f"Removed temporary file {gz_filename}")
        print("UniProt Swiss-Prot database download and extraction completed successfully!")
        return True
        
    except Exception as e:
        print(f"Error occurred: {str(e)}")
        if os.path.exists(gz_filename):
            os.remove(gz_filename)
        if os.path.exists(output_file):
            os.remove(output_file)
        return False


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Download UniProt protein database")
    parser.add_argument("-f", "--force", action="store_true", 
                       help="Force re-download even if file exists")
    parser.add_argument("-o", "--output", default="uniprots.pep",
                       help="Output filename (default: uniprots.pep)")
    
    args = parser.parse_args()
    download_uniprot_pep(force_download=args.force, output_file=args.output)


if __name__ == "__main__":
    main()

