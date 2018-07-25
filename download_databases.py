from __future__ import print_function
import os
import sys
import yaml
import argparse

def create_dir_if_absent(path):
    if not os.path.isdir(path):
        try:
            os.system("mkdir " + path)
            print("LOG: created directory: " + path, file=sys.stderr)
        except Exception as err_msg:
            sys.exit("ERROR: cannot create {}: {}".format(path, err_msg))

def download_database(db_dir, db_name):
    db_url = databases["url_" + db_name]
    db_subdir = os.path.join(db_dir, db_name)
    create_dir_if_absent(db_subdir)
    for db_file in databases['files_'+db_name]:
        db_file_url = db_url + '/' + db_file
        db_file_path = os.path.join(db_subdir, db_file)
        if (os.path.isfile(db_file_path) or 
            os.path.isfile(db_file_path.replace('.gz', ''))): # already exists
            print("LOG: {} already exists.".format(db_file_path), file=sys.stderr)
            continue
        else: # files absent
            cmd = "wget -P {db_subdir} {db_file_url}".format(
                db_subdir=db_subdir, db_file_url=db_file_url)
            try:
                os.system(cmd)
            except Exception as err_msg:
                sys.exit("ERROR: cannot fully download {}: {}".format(
                    db_file_path, err_msg))

def process_database(db_dir, db_name, args):
    db_subdir = os.path.join(db_dir, db_name)
    if not os.path.exists(os.path.join(db_subdir, "COMPLETE")):
        if db_name == "EVS":
            bash_file = os.path.join(db_dir, "process_{}.sh".format(db_name))
            os.system('bash {}'.format(bash_file))
        else:
            smk_file = os.path.join(db_dir, "process_{}.smk".format(db_name))
            cmd = ("snakemake --config configfile={} "
                "-s {}".format(args.file, smk_file))
            os.system(cmd)
    else:
        print("LOG: Database for {} already processed.".format(db_name))

if __name__ == "__main__":
    script_name = "rainbow"
    parser = argparse.ArgumentParser(description="Download databases "
        "required to run {}".format(script_name))
    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('-f', '--file', required=True, 
        help='input database-to-download yaml')
    args = parser.parse_args()
    with open(args.file, "r") as stream:
        try:
            databases = yaml.load(stream)
        except yaml.YAMLError as err_msg:
            sys.exit("ERROR: cannot load {}: {}".format(args.file, err_msg))

    db_dir = os.path.abspath(databases['database_dir'])
    create_dir_if_absent(db_dir)
    db_names = sorted([x.replace('url_', '') 
        for x in databases.keys() if x.startswith('url_')])
    for db_name in db_names:
        download_database(db_dir, db_name)
        process_database(db_dir, db_name, args)

