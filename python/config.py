import os

class Config:
    def __init__(self, cpus=1, distinct_count=1, min_len=120,
                    db_path="", db_name=""):
        self.cpus = cpus
        self.distinct_count = distinct_count
        self.min_len = min_len
        self.db_path = db_path
        self.db_name = db_name

    def print(self):
        print(self.cpus, self.distinct_count, self.min_len, self.db_name)

    def format(self):
        db = os.path.join(self.db_path, self.db_name)
        params = f"-p {self.cpus} -k {self.distinct_count} --min-hitlen {self.min_len} -x {db}"
        return params
