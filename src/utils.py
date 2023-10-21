from typing import List




def make_file_name(
	components:List[str]
	) -> str:
	# assume components[0] is fname, components[-1]
	base_dir = components[0]
	base_dir = base_dir if base_dir[-1]!="/" else base_dir[:-1]
	print(base_dir)