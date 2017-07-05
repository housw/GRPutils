#!/bin/bash


GRP_FILES=./grp_files


# --------------------
# summarize grp files
# --------------------

summarize_file=summarize_grp_files.txt

:<< 'COMMENT'
if [ ! -e "$summarize_file" ] ; then
    touch "$summarize_file"
else
    echo "$summarize_file" exists, please remove/rename it first !!
    exit
fi

for f in `ls $GRP_FILES`; do
    echo "summarizing file: $f" >> $summarize_file
    ./grptools summarize $GRP_FILES/$f >> $summarize_file
    echo -e "### ----------------\n" >> $summarize_file
done
COMMENT



# --------------------
# normalize grp files
# --------------------
# for each condition, normalize the 
# dRNA-Seq and minus libraries to the largest dRNAd-Seq 
# library size (fwd_tss+rev_tss)

NORM_FILES="01_normalize"
mkdir -p $NORM_FILES

:<< 'COMMENT'
for f in `ls $GRP_FILES`; do

    # normalize chromosome
    if [[ $f =~ "_TE101_" ]]; then
        case "$f" in 
            43ppt_*)
                echo "normalizing chromosome grp at 43ppt: $f"
                ./grptools normalize $GRP_FILES/$f 1984283 ;;
            37ppt_*)
                echo "normalizing chromosome grp at 37ppt: $f"
                ./grptools normalize $GRP_FILES/$f 1145614 ;;
            *)
                echo "$f wouldn't be normalized!";;
        esac

    # normalize plasmid
    elif [[ $f =~ "_pTE101_" ]]; then
        case "$f" in 
            43ppt_*)
                echo "normalizing plasmid grp at 43ppt: $f"
                ./grptools normalize $GRP_FILES/$f 11375 ;;
            37ppt_*)
                echo "normalizing plasmid grp at 37ppt: $f"
                ./grptools normalize $GRP_FILES/$f 8033 ;;
            *)
                echo "$f wouldn't be normalized!";;
        esac
    else
        echo "not included"
    fi
done
mv *_normalized.grp $NORM_FILES
COMMENT


# --------------------
# subtract grp files
# --------------------

SUBT_FILES="02_subtract"
mkdir -p $SUBT_FILES

:<< 'COMMENT'
for f in `ls $NORM_FILES`; do

    if [[ $f =~ "_dRNA_" ]]; then
        echo $f
        dRNA=$NORM_FILES/$f
        minus=${dRNA/_dRNA_/_minus_}
        ./grptools subtract $dRNA $minus -f sqrt  # sqrt transform second array
    fi
done
mv *_normalized_subtracted.grp $SUBT_FILES
COMMENT


# --------------------
# intersect grp files
# --------------------

INTS_FILES="03_intersect"
mkdir -p $INTS_FILES

:<< 'COMMENT'
for f in `ls $SUBT_FILES`; do

    if [[ $f =~ "_S5_" ]]; then
        echo $f
        S5=$SUBT_FILES/$f
        S6=${S5/_S5_/_S6_}
        ./grptools intersect $S5 $S6
    elif [[ $f =~ "_S3_" ]]; then
        echo $f
        S3=$SUBT_FILES/$f
        S4=${S3/_S3_/_S4_}
        ./grptools intersect $S3 $S4
    else
        echo $f
    fi
done
mv *_intersected.grp $INTS_FILES
if [ $? -ne 0 ]; then
    echo -e "\tFiles not exist, did you comment the commands out?!"
fi
COMMENT



# -----------------------------
# merge grp files by sumup mode
# -----------------------------
SUMUP_FILES="04_sumup"
mkdir -p $SUMUP_FILES

:<< 'COMMENT'
for f in `ls $INTS_FILES`; do

    if [[ $f =~ "_S5_" ]]; then
        echo $f
        S5=$INTS_FILES/$f
        S6=${S5/_S5_/_S6_}
        ./grptools merge -m sumup $S5 $S6
    elif [[ $f =~ "_S3_" ]]; then
        echo $f
        S3=$INTS_FILES/$f
        S4=${S3/_S3_/_S4_}
        ./grptools merge -m sumup $S3 $S4
    else
        echo $f
    fi
done
mv *_sumup.grp $SUMUP_FILES
COMMENT


# --------------------------------
# aggregate tss of sumup grp files
# --------------------------------
AGG_FILES="05_aggregated"
mkdir -p $AGG_FILES

:<< 'COMMENT'
for f in `ls $SUMUP_FILES`; do
    echo $f
    ./grptools aggregate $SUMUP_FILES/$f
done
mv *_aggregated.grp $AGG_FILES
COMMENT


# -----------------------------------
# compareTSS to merge 43ppt and 37ppt
# -----------------------------------
INTEGRATE_FILES="06_integrate_TSS"
mkdir -p $INTEGRATE_FILES

:<< 'COMMENT'

for replicon in _TE101_ _pTE101_ ; do
    echo $replicon
    ./grptools compareTSS `ls $AGG_FILES/*${replicon}*` -t
done
mv *_compareTSS_*.grp $INTEGRATE_FILES

COMMENT


# ---------------------------
# get tssTable from grp files
# ---------------------------
TAB_FILES='07_TSS_Tables'
mkdir -p $TAB_FILES

:<< 'COMMENT'
for f in `ls $INTEGRATE_FILES`; do
    echo $f
    ./tptools grp2TssTable $INTEGRATE_FILES/$f
done
mv *_grp2tss.tab $TAB_FILES
COMMENT

# ------------
# classify TSS
# ------------

TSS_FILES='08_classified_TSS'
mkdir -p $TSS_FILES

:<< 'COMMENT'
## 43 ppt and 37 ppt combined
./tptools classify ./$TAB_FILES/37ppt_dRNA_fwd_S3_cutadapter_t20l20_uniClustered_nonrRNA_us_mapped_to_TE101_id95_normalized_subtracted_intersected_sumup_aggregated_compareTSS_max_grp2tss.tab \
                   ../genomes/TE101.gff \
                   ./$INTEGRATE_FILES/37ppt_dRNA_fwd_S3_cutadapter_t20l20_uniClustered_nonrRNA_us_mapped_to_TE101_id95_normalized_subtracted_intersected_sumup_aggregated_compareTSS_max.grp \
                   -r /home/hou/AlteroTSS/segemehl_sam_grp/grp_files/Mischung_rdm_cutadapter_t20l20_uniClustered_nonrRNA_us_mapped_to_TE101_id95.grp

./tptools classify ./$TAB_FILES/37ppt_dRNA_fwd_S3_cutadapter_t20l20_uniClustered_nonrRNA_us_mapped_to_pTE101_id95_normalized_subtracted_intersected_sumup_aggregated_compareTSS_max_grp2tss.tab \
                   ../genomes/pTE101a.gff \
                   ./$INTEGRATE_FILES/37ppt_dRNA_fwd_S3_cutadapter_t20l20_uniClustered_nonrRNA_us_mapped_to_pTE101_id95_normalized_subtracted_intersected_sumup_aggregated_compareTSS_max.grp \
                   -r /home/hou/AlteroTSS/segemehl_sam_grp/grp_files/Mischung_rdm_cutadapter_t20l20_uniClustered_nonrRNA_us_mapped_to_pTE101_id95.grp                   
mv *_classified.tss $TSS_FILES
COMMENT


# ----------
# filter TSS
# ----------
FILTER_TSS='09_filtered_TSS'
mkdir -p $FILTER_TSS

:<< 'COMMENT'
for f in `ls $TSS_FILES`; do
    echo $f
    if [[ -d $TSS_FILES/$f ]]; then
        echo "$f is a directory, do nothing"
    else
        echo "filtering $f ..." 
        ./tptools filterTss $TSS_FILES/$f -e 10 -r 0.5 -c 0.5 -t 0.3
    fi
done
mv *_filtered_*.tss $FILTER_TSS
COMMENT


# ----------------------------
# aggregate raw dRNA grp files
# ----------------------------
AGG_RAW_FILES="10_aggregate_raw"
mkdir -p $AGG_RAW_FILES

:<< 'COMMENT'
for f in `ls $GRP_FILES`; do
    if [[ $f =~ "_dRNA_" && $f =~ "_TE101_" ]]; then
        echo $f
        ./grptools aggregate $GRP_FILES/$f -g $INTEGRATE_FILES/37ppt_dRNA_fwd_S3_cutadapter_t20l20_uniClustered_nonrRNA_us_mapped_to_TE101_id95_normalized_subtracted_intersected_sumup_aggregated_compareTSS_max.grp
    elif [[ $f =~ "_dRNA_" && $f =~ "_pTE101_" ]]; then
        echo $f
        ./grptools aggregate $GRP_FILES/$f -g $INTEGRATE_FILES/37ppt_dRNA_fwd_S3_cutadapter_t20l20_uniClustered_nonrRNA_us_mapped_to_pTE101_id95_normalized_subtracted_intersected_sumup_aggregated_compareTSS_max.grp
    fi
done
mv *_aggregated.grp ./$AGG_RAW_FILES
COMMENT


# --------------
# map grp to tss
# --------------
MAP_TSS='11_MAP_TSS'
mkdir -p $MAP_TSS

:<< 'COMMENT'
./tptools mapGrp2Tss `ls $AGG_RAW_FILES/*_TE101_*_aggregated.grp` $FILTER_TSS/37ppt_dRNA_fwd_S3_cutadapter_t20l20_uniClustered_nonrRNA_us_mapped_to_TE101_id95_normalized_subtracted_intersected_sumup_aggregated_compareTSS_max_grp2tss_classified_filtered_10.00.tss

./tptools mapGrp2Tss `ls $AGG_RAW_FILES/*_pTE101_*_aggregated.grp` $FILTER_TSS/37ppt_dRNA_fwd_S3_cutadapter_t20l20_uniClustered_nonrRNA_us_mapped_to_pTE101_id95_normalized_subtracted_intersected_sumup_aggregated_compareTSS_max_grp2tss_classified_filtered_10.00.tss
mv *_mapped.tab $MAP_TSS
COMMENT
