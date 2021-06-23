import os
from os import path
from itertools import islice 

def parse_dir(input_directories, forme  = ".tsv"):

	path_list = sorted([d for d in input_directories.iterdir() if d.is_dir()])
	paths = []

	for i in range(len(path_list)):
		input_directorie = path_list[i]
		path = sorted(list(input_directorie.glob('**/*.' + forme)))
		
		for j in range(len(path)):
			paths.append(path[j])

	return(path_list, paths)