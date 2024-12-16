import requests
import tempfile
from zipfile import ZipFile


DB_URL = "https://figshare.com/ndownloader/articles/19679595/versions/2"


def init_setup_db(output_dir="default_db"):

    data = requests.get(DB_URL)
    # with open(temp_file, 'wb') as file:
    temp_file = tempfile.TemporaryFile()
    temp_file.write(data.content)

    with ZipFile(temp_file, 'r') as zipObj:
        zipObj.extractall(output_dir)

    return True
