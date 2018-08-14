## Install
Download required softwares and databases for running RainBow.

    git clone https://github.com/rg-korea/rainbow.git
    cd rainbow
    bash setup.sh

## Activate the installed conda environment
The install file created a conda environment called rainbow.
You have to activate this environment before using rainbow.

    source activate rainbow

## Copy and edit the config .ini file
Edit the config file so that the parameters and directory names fit your 
analysis needs. The recommended method is that you copy the example.ini file 
from the `example_data` folder, and edit the contents, as follows.

    cp example_data/example.ini my_config.ini

Note that when editing the `data_dir` of the config file, you should make sure 
that all fastq files are in fastq(.gz) format, and the files should end in either
`_1.fastq.gz`, `_1.fastq`, `_2.fastq.gz` or `_2.fastq`.

## Run
After getting your own copy of the config file, now you can run the 
rainbow pipeline.

    python rainbow.py -i /path/to/your/my_config.ini

## Third party software/database license
RainBow contains files not developed by the Rare Genomics Institute Korea.
Files in the subdirectory bin/ and files downloaded by setup.py follow their
own licenses. We recommend you read them, since not all of them follow the same 
license as RainBow.

## License
RainBow follows The MIT License (MIT). This license applies to all parts of this
software except the third party software/databases documented in 
THIRD-PARTY-NOTICES.txt. 
For further information, we recommend you read our LICENSE.txt or the contents 
in https://opensource.org/licenses/MIT

