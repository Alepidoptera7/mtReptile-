line_list = []
with open("percent_ids.txt") as in_file:
    for line in in_file:
        stripped_line = line.replace("(","").replace(")","")
        split_line = stripped_line.split(",")
        line_list.append(split_line)

percent_id_dict = {line[0][:15]: [] for line in line_list}

for line in line_list:
    if line[0][:15] in percent_id_dict.keys():
        percent_id_dict[line[0][:15]].append((line[0].split("a")[1].split("-")[0], line[1].strip("\n")))


for key in percent_id_dict.keys():

    percent_id1 = float(percent_id_dict[key][0][1])
    percent_id2 = float(percent_id_dict[key][1][1])

    if percent_id1 > percent_id2:
        print(key, ": ", percent_id_dict[key][0])

    else: #percent_id_dict[key][0][1] < percent_id_dict[key][1][1]:
        print(key, ": ", percent_id_dict[key][1])

    if percent_id_dict[key][0][1] == percent_id_dict[key][1][1]:
        print(key, ": ", percent_id_dict[key])