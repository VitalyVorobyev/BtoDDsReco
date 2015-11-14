#! /bin/csh -f

setenv USE_GRAND_REPROCESS_DATA 1
source /sw/belle/local/etc/cshrc_general

basf <<EOF >& logs/log_sigmc_dds_43-0.txt
path create main
path add_module main fix_mdst
path create analysis
path add_condition main <=:0:KILL
path add_condition main >:0:analysis
path add_module analysis b2dds
#path add_module analysis user_index
path add_condition analysis <=:0:KILL
initialize

nprocess set 0
#table save mdst_all
#table save evtvtx_all
#table save hepevt_all
#table save belletdf_all
##table save evtcls_hadron_info
##table save evtcld_hadronic_flag
#table save reccdc_timing

module put_parameter b2dds ofile\"/gpfs/home/belle/vitaly/work/BtoDDs/tuples/b2dds_sigmc_43-0.root"
module put_parameter b2dds mode\1
# 0 -> Data, 1 -> Signal MC, 2 -> Generic MC
module put_parameter b2dds ntuple_flag\1

process_event /gpfs/home/belle/vitaly/work/mixing/MCProd/gsim/mdst/BtoDDs/evtgen_exp_43_BtoDDs-0.mdst

EOF
