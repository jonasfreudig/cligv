import os
import sys

def clean_null_bytes(filename):
    with open(filename, 'rb') as f:
        content = f.read()
    
    # Remove null bytes
    clean_content = content.replace(b'\x00', b'')
    
    cleaned_filename = f"clean_{filename}"
    with open(cleaned_filename, 'wb') as f:
        f.write(clean_content)
    
    print(f"Cleaned {filename}")
    return cleaned_filename

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python clean.py file1.py file2.py ...")
        sys.exit(1)
    
    for filename in sys.argv[1:]:
        cleaned_file = clean_null_bytes(filename)
        os.rename(cleaned_file, filename)
        print(f"Replaced {filename}")
