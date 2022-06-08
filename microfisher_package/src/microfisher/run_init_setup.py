import sys
import requests
import tempfile
from zipfile import ZipFile


def init_setup_db(output_dir = "default_db"):

    URL = "https://figshare.com/ndownloader/articles/19679595/versions/1"
    data = requests.get(URL)

    # with open(temp_file, 'wb') as file:
    temp_file = tempfile.TemporaryFile()
    temp_file.write(data.content)

    with ZipFile(temp_file, 'r') as zipObj:
       zipObj.extractall(output_dir)

    return True

#
# if __name__ == "__main__":
#     main()
