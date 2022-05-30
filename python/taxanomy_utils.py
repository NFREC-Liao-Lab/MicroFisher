from taxanomy import get_desired_taxa_ranks

def summarise_data(parsed_data, rank):
    results = dict()
    for info in parsed_data:
        try:
            key = info[rank]
        except KeyError as e:
            continue
        count = int(info["numReads"])
        try:
            results[key].append(count)
        except KeyError as e:
            results[key] = [count]
    summary = {k:sum(v) for k, v in results.items()}
    return(summary)


def process_report(all_lines, heading):
    parsed_data = list()
    for t in all_lines:
        ts = t.strip().split("\t")
        info = dict(zip(heading, ts))
        rank_data = get_desired_taxa_ranks(info["taxID"])
        info.update(rank_data)
        parsed_data.append(info)
    return(parsed_data)



def normalise_data(data_each):
    total = sum(data_each.values())
    data_n = {k:(v/total) for k, v in data_each.items()}
    return(data_n)


#
# report_files = "example.report.tsv"
# report_files = ["eg1.report.tsv", "eg2.report.tsv"]

def parse_report_files(report_files):
    parsed_data = dict()
    for f in report_files:
        with open(f, "r") as report:
            heading = report.readline().strip().split("\t")
            all_lines = report.readlines()
            parsed_data[f] = process_report(all_lines, heading)
    return(parsed_data)
