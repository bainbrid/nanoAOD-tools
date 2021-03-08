function run_test()
{
    cd ~
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    export SCRAM_ARCH=slc6_amd64_gcc700 || return 1
    scramv1 project CMSSW CMSSW_10_2_18 || return 1
    cd CMSSW_10_2_18/src || return 1
    eval `scram runtime -sh` || return 1
    rsync -r /scripts/* LLPReco  || return 1
    scram b || return 1
    mkdir -p PhysicsTools/NanoAODTools
    rsync -r --stats /scripts/ PhysicsTools/NanoAODTools/. || return 1
    scram b || return 1
    echo "--- Test HNL script ---"
    # add data test
    python PhysicsTools/NanoAODTools/processors/HNL.py --year 2016 --leptons 1 --input=https://github.com/LLPDNNX/test-files/raw/master/nanoaod/Moriond17_aug2018_miniAODv3_HNL_nanoAODv3.root . || return 1
    python PhysicsTools/NanoAODTools/processors/HNL.py --year 2016 --leptons 2 --input=https://github.com/LLPDNNX/test-files/raw/master/nanoaod/Moriond17_aug2018_miniAODv3_HNL_nanoAODv3.root . || return 1
    python PhysicsTools/NanoAODTools/processors/HNL.py --year 2016 --leptons 1 --isData --input=https://github.com/LLPDNNX/test-files/raw/master/nanoaod/Moriond17_aug2018_miniAODv3_data_nanoAODv3.root . || return 1
    python PhysicsTools/NanoAODTools/processors/HNL.py --year 2016 --leptons 2 --isData --input=https://github.com/LLPDNNX/test-files/raw/master/nanoaod/Moriond17_aug2018_miniAODv3_data_nanoAODv3.root . || return 1
    python PhysicsTools/NanoAODTools/processors/HNL.py --year 2016 --leptons 1 --notagger --input=https://github.com/LLPDNNX/test-files/raw/master/nanoaod/Moriond17_aug2018_miniAODv3_data_nanoAODv3.root . || return 1
    python PhysicsTools/NanoAODTools/processors/HNL.py --year 2016 --leptons 2 --nobdt --input=https://github.com/LLPDNNX/test-files/raw/master/nanoaod/Moriond17_aug2018_miniAODv3_data_nanoAODv3.root . || return 1
    python PhysicsTools/NanoAODTools/processors/HNL.py --year 2016 --leptons 2 --notagger --notrigger --input=https://github.com/LLPDNNX/test-files/raw/master/nanoaod/Moriond17_aug2018_miniAODv3_data_nanoAODv3.root . || return 1
    python PhysicsTools/NanoAODTools/processors/HNL.py --year 2016 --leptons 2 --notagger --noiso --input=https://github.com/LLPDNNX/test-files/raw/master/nanoaod/Moriond17_aug2018_miniAODv3_data_nanoAODv3.root . || return 1
}

run_test
