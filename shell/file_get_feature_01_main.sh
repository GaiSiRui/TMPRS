>$time_point.log
>${threshold}.log
num=0
logistic="False"
clumped="False"
calculate_mod="False"
filed="False"
testmod="False"
time_point=`date +%Y%m%d%H%M%S`
prsice="False"
group="False"
test_group="False"
conv="False"
validationfile="validationfile"
outpercent=1
validationpercent=1
threshold_t=0.05
begin_site=10
end_site=200
step=10
minsite=0.00000005
takestep=0.00005
symbol=0
testfile="testfile"
merge="False"
summary="False"
validationmod="False"
with_origin="False"

source ./validfloat.sh
source ./same_fam.sh
source ./relu.sh

ARGS=`getopt -o hi:v:B:E:s:bkc:t:H:gGVI:L:T:K:mSw --long help,input-file:,validation-file:,begin-site:,end-site:,site-step:,beta,keep-t-value,calculate-mod:,test-file:,threshold-t:,group-input,group-test,group-validation,interval:,lower:,thread:,ethnic_w.txt,merge,summary,with_origin -- "$@"`
if [ $? != 0 ];then
        echo "Termination..." >&2
        exit 1
fi

eval set -- "${ARGS}"
while true
do
        case "$1" in
                -h|--help)
                        echo -e "-i|--input-file : The name of input_file, without suffix\n\n-b|--beta : Base association file\n\n-c|--calculate-mod : The calculate-mods. The chosenable label include: \n\n\t\t\t\tOrigin\n\n\t\t\t\tMean-value\n\n\t\t\t\tT-compare\n\n\t\t\t\tBeta-compare\n\n-g|--group-input : Input list mode\n\n\t\t\t\tHere, please make sure all the input file have the same .fam file.\n\n-G|--group-test : Test list mode\n\n\t\t\t\tHere, please make sure all the test file have the same .fam file.\n\n-V|--group-validation : Validation list mode\n\n\t\t\t\tHere, please make sure all the validation file have the same .fam file.\n\n-v|--validation-file : The name of validation file, if did not use this label, the system will use the last ten percent of the input dataset as the validation set\n\n-t|--test-file : The name of test file, if did not use this label, the system will use the last ten percent of the input dataset as the test set. If both of this label and the \"-v\" label have not used, the system will use the last ten percent of the input dataset as the test set, and next last ten percent of the input dataset as the validation set\n\n-B|--begin-site : The lower threshold of the split position of the system. The default value is 0. The value is no less than 0. Notice that the number can set a little higher when useing the t-split mode.\n\n-E|--end-site : The higher threshold of the split position of the system. The default value is 200. Notice that the number can set a little lower when useing the t-split mode.\n\n-s|--site-step : The threshold's step, too small value may cause the code's running step slow down. The value should no lower than (end_site - begin_site). The dafault value is 10\n\n-I|--interval : The step size of the threshold. The default is 5e-05\n\n-L|--lower : The starting p-value threshold. The default is 5e-08\n\n-h|--help : Print this help\n\n-T|--Thread : The number of thread, recommended no larger than the number of (end_site - begin_site) / site_step\n\n-H|--threshold_t : The highest p value of the snps included into the system. The default is 0.05\n\n-K|--keep_list: The list of people keep in all the process. The file need two column, the first one is family ID and the second is individual ID.\n\n-m|--merge : Merge the result of additive and non-additive calculate model\n\n-w|--with_origin : Calculate origin model in the same time calculating the other model."
                        exit 0
                        ;;

                -i|--input-file)
                        echo "input-file\t"$2"\n"
                        inputfile=$2
                        keep_file=$2
                        shift 2
                        ;;

                -v|--validation-file)
                        echo "V\t"$2"\n"
                        if [[ $2 =~ "\." ]]; then
                                echo "The validation-file cannot include the character \".\""
                                exit 1
                        fi
                        validationmod="validation"
                        validationfile=$2
                        shift 2
                        ;;

                -B|--begin-site)
                        if [[ `echo"$2<0"|bc` -eq 1 ]]; then
                                echo "The begin-site cannot smaller the 0 \".\""
                                exit 1
                        elif [ validfloat $2 ] 2>/dev/null; then
                                echo "The begin-site can only be a number."
                                exit 1
                        fi
                        echo "begin-site\t"$2"\n"
                        begin_site=$2
                        shift 2
                        ;;

                -E|--end-site)
                        if [[ `echo"$2<$begin_site"|bc` -eq 1 ]]; then
                                echo "The end-site cannot smaller than the begin-site \".\""
                                exit 1
                        elif [[ `echo"$2=$begin_site"|bc` -eq 1 ]]; then
                                echo "The end-site cannot equal to the begin-site \".\""
                                exit 1
                        elif [ validfloat $2 ] 2>/dev/null; then
                                echo "The end-site can only be a number."
                                exit 1
                        fi
                        echo "end-site\t"$2"\n"
                        end_site=$2
                        shift 2
                        ;;

                -s|--site-step)
                        if [ validfloat $2 ] 2>/dev/null; then
                                echo "The site-step can only be a number."
                                exit 1
                        fi
                        echo "s\t"$2"\n"
                        step=$2
                        shift 2
                        ;;

                -b|--beta)
                        echo "--beta\n"
                        logistic="True"
                        shift
                        ;;

                -k|--keep-t-value)
                        echo "keep-t-value\n"
                        keep_t_value="True"
                        shift
                        ;;

                -c|--calculate-mod)
                        calculatetype $2

                        if [ $? -eq 1 ] 2>/dev/null; then
                                echo "The site-step can only be one of the belows: \"False\", \"Origin\", \"Compare\", \"True\"."
                                exit 1
                        else
                                echo "OK"
                        fi

                        echo "calculate-mod:\t"$2"\n"
                        calculate_mod=$2
                        shift 2
                        ;;

                -t|--test-file)
                        echo "V\t"$2"\n"
                        if [[ $2 =~ "\." ]]; then
                                echo "The test-file cannot include the character \".\""
                                exit 1
                        fi
                        echo "test-file\t"$2"\n"
                        testmod="test"
                        testfile=$2
                        shift 2
                        ;;

                -t|--side-file)
                        echo "V\t"$2"\n"
                        if [[ $2 =~ "\." ]]; then
                                echo "The side-file cannot include the character \".\""
                                exit 1
                        fi
                        echo "side-file\t"$2"\n"
                        sidefile=$2
                        shift 2
                        ;;

                -H|--threshold-t)
                        if [[ `echo"$2<0"|bc` -eq 1 ]]; then
                                echo "The end-site cannot smaller than 0 \".\""
                                exit 1
                        elif [[ `echo"$2>1"|bc` -eq 1 ]]; then
                                echo "The end-site cannot bigger than 1 \".\""
                                exit 1
                        elif [ validfloat $2 ] 2>/dev/null; then
                                echo "The end-site can only be a number."
                                exit 1
                        fi
                        echo "threshold-t\t"$2"\n"
                        threshold_t=$2
                        shift 2
                        ;;

                -g|--group-input)
                        echo "group-input\n"
                        if [ same_fam $2 ] 2>/dev/null; then
                                echo "Please the list of input file have the same fam file, the system cannot handle the input file with differ fam file, please try to merge the bfile by plink first, sorry."
                                exit 1
                        fi
                        group="True"
                        shift
                        ;;

                -G|--group-test)
                        echo "group-test\n"
                        if [ same_fam $2 ] 2>/dev/null; then
                                echo "Please the list of test file have the same fam file, the system cannot handle the test file with differ fam file, please try to merge the bfile by plink first, sorry."
                                exit 1
                        fi
                        test_group="True"
                        shift
                        ;;

                -V|--group-validation)
                        echo "group-validation\n"
                        if [ same_fam $2 ] 2>/dev/null; then
                                echo "Please the list of validation file have the same fam file, the system cannot handle the validation file with differ fam file, please try to merge the bfile by plink first, sorry."
                                exit 1
                        fi
                        validation_group="True"
                        shift
                        ;;

                -I|--interval)
                        if [[ `echo "$2<0"|bc` -eq 1 ]]; then
                                echo "The interval cannot smaller than 0 \".\""
                                exit 1
                        elif [[ `echo "$2>1"|bc` -eq 1 ]]; then
                                echo "The interval cannot bigger than 1 \".\""
                                exit 1
                        elif [ validfloat $2 ] 2>/dev/null; then
                                echo "The interval can only be a number."
                                exit 1
                        fi
                        echo "interval\t"$2"\n"
                        takestep=$2
                        shift 2
                        ;;

                -L|--lower)
                        if [[ `echo "$2<0"|bc` -eq 1 ]]; then
                                echo "The lower cannot smaller than 0 \".\""
                                exit 1
                        elif [[ `echo "$2>1"|bc` -eq 1 ]]; then
                                echo "The lower cannot bigger than 1 \".\""
                                exit 1
                        elif [ validfloat $2 ] 2>/dev/null; then
                                echo "The lower can only be a number."
                                exit 1
                        fi
                        echo "lower\t"$2"\n"
                        minsite=$2
                        shift 2
                        ;;

                -T|--thread)
                        if [ validint $2 ] 2>/dev/null; then
                                echo "The thread can only be a integer."
                                exit 1
                        fi
                        echo "thread\t"$2"\n"
                        thread=$2
                        shift 2
                        ;;

                -m|--merge)
                        echo "merge\n"
                        merge="True"
                        shift
                        ;;

                -S|--summary)
                        echo "summary\n"
                        summary="True"
                        shift
                        ;;
                        
                -w|--with_origin)
                        echo "with_origin\n"
                        with_origin="True"
                        shift
                        ;;

                -K|--keep_list)
                        echo "K\t"$2"\n"
                        if [[ $2 =~ "\." ]]; then
                                echo "The keep-file cannot include the character \".\""
                                exit 1
                        fi
                        echo "keep-file\t"$2"\n"
                        keep_file=$2
                        shift 2
                        ;;

                --)
                        shift
                        break
                        ;;

                *)
                        echo "Wrong Options"
                        exit 7;;
        esac
done

echo "285"
if [[ $summary == "False" ]];then
        >${testfile}.summary_and_bim
        if [[ $testmod == "test" ]] && [[ $validationmod == "validation" ]];then
                fam_num=$(echo $input | awk 'END{print NR}' $inputfile.fam)
                threshold_t=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c "import t_to_p;print(t_to_p.t_to_p("${threshold_t},${fam_num}"))"`
                echo "threshold_t\t"$threshold_t
        elif [[ $testmod != "test" ]] && [[ $validationmod == "validation" ]];then
                fam_num=$(echo $input | awk 'END{print NR*9/10}' $inputfile.fam)
                threshold_t=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c "import t_to_p;print(t_to_p.t_to_p("${threshold_t},${fam_num}"))"`
                echo "threshold_t\t"$threshold_t
                echo "296"
        elif [[ $testmod == "test" ]] && [[ $validationmod != "validation" ]];then
                fam_num=$(echo $input | awk 'END{print NR*9/10}' $inputfile.fam)
                threshold_t=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c "import t_to_p;print(t_to_p.t_to_p("${threshold_t},${fam_num}"))"`
                echo "threshold_t\t"$threshold_t
        else
                fam_num=$(echo $input | awk 'END{print NR*4/5}' $inputfile.fam)
                fam_num_validation=$(echo $input | awk 'END{print NR*9/10}' $inputfile.fam)
                threshold_t=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c "import t_to_p;print(t_to_p.t_to_p("${threshold_t},${fam_num}"))"`
                echo "threshold_t\t"$threshold_t
        fi
        echo "307"
        
        if [ $test_group = "True" ] ;then
                testfile_num=$(echo $input | awk 'END{print NR}' ${testfile})
        fi
        if [ $group = "True" ] ;then
                inputfile_num=$(echo $input | awk 'END{print NR}' ${inputfile})
                inputfile_site=$(echo $input | awk 'NR==1{print $0}' ${inputfile})
                awk 'NR>1{print $0"_few"}' ${inputfile} > ${inputfile}_few
                >${testfile}.list
                >${inputfile}.list
                echo "318"
                >${validationfile}.list
                if [[ $testmod != "test" ]] && [[ $validationmod == "validation" ]];then
                        awk -v line_a=$fam_num 'NR>line_a&&NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${inputfile_site}.fam ${keep_file} >> ${testfile}.list
                        awk -v line_a=$fam_num 'NR>1&&NR<=line_a{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${inputfile_site}.fam ${keep_file} >> ${inputfile}.list
                        awk 'NR>1&&NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${validationfile}.fam ${keep_file} >> ${validationfile}.list
                elif [[ $testmod == "test" ]] && [[ $validationmod != "validation" ]];then
                        awk -v line_a=$fam_num 'NR>line_a&&NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${inputfile_site}.fam ${keep_file} >> ${validationfile}.list
                        awk -v line_a=$fam_num 'NR>1&&NR<=line_a&&NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${inputfile_site}.fam ${keep_file} >> ${inputfile}.list
                        awk 'NR>1&&NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${testfile}.fam ${keep_file} >> ${testfile}.list
                elif [[ $testmod == "test" ]] && [[ $validationmod == "validation" ]];then
                        echo "329"
                        awk 'NR>1&&NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${validationfile}.fam ${keep_file} >> ${validationfile}.list
                        awk 'NR>1&&NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${inputfile_site}.fam ${keep_file} >> ${inputfile}.list
                        awk 'NR>1&&NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${testfile}.fam ${keep_file} >> ${testfile}.list
                else
                        awk -v line_a=$fam_num_validation 'NR>line_a&&NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${inputfile_site}.fam ${keep_file} >> ${testfile}.list
                        awk -v line_a=$fam_num 'NR>1&&NR<=line_a&&NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${inputfile_site}.fam ${keep_file} >> ${inputfile}.list
                        awk -v line_a=$fam_num -v line_b=$fam_num_validation 'NR>line_a&&NR<=line_b&&NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${inputfile_site}.fam ${keep_file} >> ${validationfile}.list
                fi
                comm -12 <(sort ${inputfile}.list|uniq) <(sort ${inputfile}.list|uniq) > time
                mv time ${inputfile}.list
                echo "340"
                comm -12 <(sort ${validationfile}.list|uniq) <(sort ${validationfile}.list|uniq) > time
                mv time ${validationfile}.list
                comm -12 <(sort ${testfile}.list|uniq) <(sort ${testfile}.list|uniq) > time
                mv time ${testfile}.list
                awk '{print $1"\t"$1"\t"$6}' ${inputfile}.fam > ${inputfile}.pheno
                awk '{print $1"\t"$1"\t"$5}' ${inputfile}.fam > ${inputfile}.sex
                >${inputfile}.snps
                >${inputfile}.assoc
                for inputfile_order in $(seq 1 $inputfile_num)
                do
                        echo "351"
                        inputfile_site=$(echo $input | awk 'NR=="'$inputfile_order'"{print $0}' ${inputfile})
                        if [ $logistic = "True" ] ;then
                                if [ $calculate_mod = "Origin" ];then
                                        plink --bfile ${inputfile_site}  --pheno ${inputfile}.pheno --keep ${inputfile}.list --geno 0.05 --maf 0.05 --out ${inputfile_site} --logistic hide-covar --covar ${sidefile} keep-pheno-on-missing-cov
                                        awk '$0~/ADD/||NR==1{print $0}' ${inputfile_site}.assoc.logistic > ${inputfile_site}.assoc.add
                                        if [ $clumped = "False" ];then
                                                plink --bfile ${inputfile_site} --keep ${inputfile}.list --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump ${inputfile_site}.assoc.logistic --clump-snp-field SNP --clump-field P
                                        fi
                                        awk 'NR==FNR{a[$3]=$5};NR>FNR{if(($2 in a)&&(a[$2]<0.5))print $2"\t"$7"\t"$8"\t"$9}' plink.clumped ${inputfile_site}.assoc.logistic > ${inputfile_site}.assoc
                                        echo "362"
                                        awk '{print $1}' ${inputfile_site}.assoc > ${inputfile_site}.snps
                                else
                                        plink --bfile ${inputfile_site}  --pheno ${inputfile}.pheno --keep ${inputfile}.list --geno 0.05 --maf 0.05 --out ${inputfile_site} --logistic dominant hide-covar --covar ${sidefile} keep-pheno-on-missing-cov
                                        awk '$0~/DOM/||NR==1{print $0}' ${inputfile_site}.assoc.logistic > ${inputfile_site}.assoc.dom
                                        plink --bfile ${inputfile_site}  --pheno ${inputfile}.pheno --keep ${inputfile}.list --geno 0.05 --maf 0.05 --out ${inputfile_site} --logistic recessive hide-covar --covar ${sidefile} keep-pheno-on-missing-cov
                                        awk '$0~/REC/||NR==1{print $0}' ${inputfile_site}.assoc.logistic > ${inputfile_site}.assoc.rec
                                        plink --bfile ${inputfile_site}  --pheno ${inputfile}.pheno --keep ${inputfile}.list --geno 0.05 --maf 0.05 --out ${inputfile_site} --logistic hide-covar --covar ${sidefile} keep-pheno-on-missing-cov
                                        awk '$0~/ADD/||NR==1{print $0}' ${inputfile_site}.assoc.logistic > ${inputfile_site}.assoc.add
                                        cat ${inputfile_site}.assoc.add ${inputfile_site}.assoc.dom ${inputfile_site}.assoc.rec > ${inputfile_site}.assoc.logistic
                                        if [ $clumped = "False" ];then
                                                echo "373"
                                                plink --bfile ${inputfile_site} --keep ${inputfile}.list --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump ${inputfile_site}.assoc.logistic --clump-snp-field SNP --clump-field P
                                        fi
                                        awk 'NR==FNR{a[$3]=$5};NR>FNR{if(($2 in a)&&(a[$2]<0.5))print $2"\t"$7"\t"$8"\t"$9}' plink.clumped ${inputfile_site}.assoc.add > ${inputfile_site}.assoc
                                        awk 'NR==FNR{a[$1]=$0};NR>FNR{if($2 in a)print a[$2]"\t"$7"\t"$8"\t"$9}' ${inputfile_site}.assoc ${inputfile_site}.assoc.dom > ${inputfile_site}.dom
                                        awk 'NR==FNR{a[$1]=$0};NR>FNR{if($2 in a)print a[$2]"\t"$7"\t"$8"\t"$9}'  ${inputfile_site}.dom ${inputfile_site}.assoc.rec > ${inputfile_site}.assoc
                                        awk '{print $1}' ${inputfile_site}.assoc > ${inputfile_site}.snps
                                fi
                        else
                                if [ $calculate_mod = "Origin" ];then
                                        plink --bfile ${inputfile_site}  --pheno ${inputfile}.pheno --keep ${inputfile}.list --geno 0.05 --maf 0.05 --out ${inputfile_site} --linear hide-covar --covar ${sidefile} keep-pheno-on-missing-cov
                                        echo "384"
                                        awk '$0~/ADD/||NR==1{print $0}' ${inputfile_site}.assoc.linear > ${inputfile_site}.assoc.add
                                        if [ $clumped = "False" ];then
                                                plink --bfile ${inputfile_site} --keep ${inputfile}.list --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump ${inputfile_site}.assoc.linear --clump-snp-field SNP --clump-field P
                                        fi
                                        awk 'NR==FNR{a[$3]=$5};NR>FNR{if(($2 in a)&&(a[$2]<0.5))print $2"\t"$7"\t"$8"\t"$9}' plink.clumped ${inputfile_site}.assoc.add > ${inputfile_site}.assoc
                                        awk '{print $1}' ${inputfile_site}.assoc > ${inputfile_site}.snps
                                else
                                        plink --bfile ${inputfile_site}  --pheno ${inputfile}.pheno --keep ${inputfile}.list --geno 0.05 --maf 0.05 --out ${inputfile_site} --linear dominant hide-covar --covar ${sidefile} keep-pheno-on-missing-cov
                                        awk '$0~/DOM/||NR==1{print $0}' ${inputfile_site}.assoc.linear > ${inputfile_site}.assoc.dom
                                        echo "395"
                                        plink --bfile ${inputfile_site}  --pheno ${inputfile}.pheno --keep ${inputfile}.list --geno 0.05 --maf 0.05 --out ${inputfile_site} --linear recessive hide-covar --covar ${sidefile} keep-pheno-on-missing-cov
                                        awk '$0~/REC/||NR==1{print $0}' ${inputfile_site}.assoc.linear > ${inputfile_site}.assoc.rec
                                        plink --bfile ${inputfile_site}  --pheno ${inputfile}.pheno --keep ${inputfile}.list --geno 0.05 --maf 0.05 --out ${inputfile_site} --linear hide-covar --covar ${sidefile} keep-pheno-on-missing-cov
                                        awk '$0~/ADD/||NR==1{print $0}' ${inputfile_site}.assoc.linear > ${inputfile_site}.assoc.add
                                        cat ${inputfile_site}.assoc.add ${inputfile_site}.assoc.dom ${inputfile_site}.assoc.rec > ${inputfile_site}.assoc.linear
                                        if [ $clumped = "False" ];then
                                                plink --bfile ${inputfile_site} --keep ${inputfile}.list --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump ${inputfile_site}.assoc.linear --clump-snp-field SNP --clump-field P
                                        fi
                                        awk 'NR==FNR{a[$3]=$5};NR>FNR{if(($2 in a)&&(a[$2]<0.5))print $2"\t"$7"\t"$8"\t"$9}' plink.clumped ${inputfile_site}.assoc.add > ${inputfile_site}.assoc
                                        awk 'NR==FNR{a[$1]=$0};NR>FNR{if($2 in a)print a[$2]"\t"$7"\t"$8"\t"$9}' ${inputfile_site}.assoc ${inputfile_site}.assoc.dom > ${inputfile_site}.dom
                                        echo "406"
                                        awk 'NR==FNR{a[$1]=$0};NR>FNR{if($2 in a)print a[$2]"\t"$7"\t"$8"\t"$9}'  ${inputfile_site}.dom ${inputfile_site}.assoc.rec > ${inputfile_site}.assoc
                                        awk '{print $1}' ${inputfile_site}.assoc > ${inputfile_site}.snps
                                fi
                        fi
                        mv plink.clumped ${inputfile_site}.clumped
                        cat ${inputfile}.clumped ${inputfile_site}.clumped > time
                        uniq time > plink.clumped
                        sort -g -k 5 plink.clumped -o time
                        sed '/^$/d' time > plink.clumped
                        cat ${inputfile}.snps ${inputfile_site}.snps > time
                        echo "417"
                        uniq time > ${inputfile}.snps
                        cat ${inputfile}.assoc ${inputfile_site}.assoc > time
                        uniq time > ${inputfile}.assoc
                done
                declare -A snp_site
                for snp_num_site in $(seq $minsite $takestep 0.05)
                do
                        snp_num=$(echo $input | awk -v a="$snp_num_site" 'NR>1{if($5>a){print NR-1;exit}}' plink.clumped)
                        snp_site[${#snp_site[@]}]=$snp_num
                        snp_num=$(echo $input | awk -v a="$snp_num_site" 'NR>1{if($5>a){print int((NR-2)/1000)*1000;exit}}' plink.clumped)
                        echo "428"
                        snp_site_down[${#snp_site_down[@]}]=$snp_num
                done
                echo ${snp_site[*]}
                echo ${snp_site_down[*]}
                input_num=$(echo $input | awk 'END{print NR}' ${inputfile}.snps)
                input_length=$((input_num/1000))
                split -l 1000 -d -a ${#input_length} ${inputfile}.snps ${inputfile}_
                split -l 1000 -d -a ${#input_length} ${inputfile}.assoc ${inputfile}.assoc_
                echo $input_length
                echo ${#input_length}
                echo "439"
                input_len=$((input_length-1))
                for inputfile_order in $(seq 1 $inputfile_num)
                do
                    inputfile_site=$(echo $input | awk 'NR=="'$inputfile_order'"{print $0}' ${inputfile})
                    plink --extract ${inputfile}.snps --bfile ${inputfile_site}  --make-bed --out ${inputfile_site}_few
                done
                inputfile_site=$(echo $input | awk 'NR==1{print $0}' ${inputfile})
                plink --noweb --bfile ${inputfile_site}_few --merge-list ${inputfile}_few --make-bed --out ${inputfile}
        else
                if [[ $testmod != "test" ]] && [[ $validationmod == "validation" ]];then
                        echo "450"
                        awk -v line_a=$fam_num 'NR>line_a{print $1"\t"$1}' ${inputfile}.fam > ${testfile}.list
                        awk -v line_a=$fam_num 'NR>1&&NR<=line_a{print $1"\t"$1}' ${inputfile}.fam > ${inputfile}.list
                        awk 'NR>1{print $1"\t"$1}' ${validationfile}.fam > ${validationfile}.list
                elif [[ $testmod == "test" ]] && [[ $validationmod != "validation" ]];then
                        awk -v line_a=$fam_num 'NR>line_a{print $1"\t"$1}' ${inputfile}.fam > ${validationfile}.list
                        awk -v line_a=$fam_num 'NR>1&&NR<=line_a{print $1"\t"$1}' ${inputfile}.fam > ${inputfile}.list
                        awk 'NR>1{print $1"\t"$1}' ${testfile}.fam > ${testfile}.list
                elif [[ $testmod == "test" ]] && [[ $validationmod == "validation" ]];then
                        awk '{print $1"\t"$1}' ${validationfile}.fam > ${validationfile}.list
                        awk '{print $1"\t"$1}' ${inputfile}.fam > ${inputfile}.list
                        awk '{print $1"\t"$1}' ${testfile}.fam > ${testfile}.list
                else
                        awk -v line_a=$fam_num_validation 'NR>line_a{print $1"\t"$1}' ${inputfile}.fam > ${testfile}.list
                        awk -v line_a=$fam_num 'NR>1&&NR<=line_a{print $1"\t"$1}' ${inputfile}.fam > ${inputfile}.list
                        awk -v line_a=$fam_num -v line_b=$fam_num_validation 'NR>line_a&&NR<=line_b{print $1"\t"$1}' ${inputfile}.fam > ${validationfile}.list
                fi
                comm -12 <(sort ${inputfile}.list|uniq) <(sort ${inputfile}.list|uniq) > time
                mv time ${inputfile}.list
                comm -12 <(sort ${validationfile}.list|uniq) <(sort ${validationfile}.list|uniq) > time
                mv time ${validationfile}.list
                comm -12 <(sort ${testfile}.list|uniq) <(sort ${testfile}.list|uniq) > time
                mv time ${testfile}.list
                echo "472"
                awk '{print $1"\t"$1"\t"$6}' ${inputfile}.fam > ${inputfile}.pheno
                awk '{print $1"\t"$1"\t"$5}' ${inputfile}.fam > ${inputfile}.sex
                if [ $logistic = "True" ] ;then
                        if [ $calculate_mod = "Origin" ];then
                                plink --bfile ${inputfile}  --pheno ${inputfile}.pheno --keep ${inputfile}.list --geno 0.05 --maf 0.05 --out ${inputfile} --logistic hide-covar --covar ${sidefile} keep-pheno-on-missing-cov
                                awk '$0~/ADD/||NR==1{print $0}' ${inputfile}.assoc.logistic > ${inputfile}.assoc.add
                                mv ${inputfile}.assoc.add ${inputfile}.assoc.logistic
                                if [ $clumped = "False" ];then
                                        echo "483"
                                        plink --bfile ${inputfile} --keep ${inputfile}.list --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump ${inputfile}.assoc.logistic --clump-snp-field SNP --clump-field P
                                fi
                                awk 'NR==FNR{a[$3]=$5};NR>FNR{if(($2 in a)&&(a[$2]<0.5))print $2"\t"$7"\t"$8"\t"$9}' plink.clumped ${inputfile}.assoc.logistic > ${inputfile}.assoc
                                awk '{print $1}' ${inputfile}.assoc > ${inputfile}.snps
                        else
                                plink --bfile ${inputfile}  --pheno ${inputfile}.pheno --keep ${inputfile}.list --geno 0.05 --maf 0.05 --out ${inputfile} --logistic dominant hide-covar --covar ${sidefile} keep-pheno-on-missing-cov
                                awk '$0~/DOM/||NR==1{print $0}' ${inputfile}.assoc.logistic > ${inputfile}.assoc.dom
                                plink --bfile ${inputfile}  --pheno ${inputfile}.pheno --keep ${inputfile}.list --geno 0.05 --maf 0.05 --out ${inputfile} --logistic recessive hide-covar --covar ${sidefile} keep-pheno-on-missing-cov
                                awk '$0~/REC/||NR==1{print $0}' ${inputfile}.assoc.logistic > ${inputfile}.assoc.rec
                                plink --bfile ${inputfile}  --pheno ${inputfile}.pheno --keep ${inputfile}.list --geno 0.05 --maf 0.05 --out ${inputfile} --logistic hide-covar --covar ${sidefile} keep-pheno-on-missing-cov
                                echo "494"
                                awk '$0~/ADD/||NR==1{print $0}' ${inputfile}.assoc.logistic > ${inputfile}.assoc.add
                                cat ${inputfile}.assoc.add ${inputfile}.assoc.dom ${inputfile}.assoc.rec > ${inputfile}.assoc.logistic
                                if [ $clumped = "False" ];then
                                        plink --bfile ${inputfile} --keep ${inputfile}.list --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump ${inputfile}.assoc.logistic --clump-snp-field SNP --clump-field P
                                fi
                                awk 'NR==FNR{a[$3]=$5};NR>FNR{if(($2 in a)&&(a[$2]<0.5))print $2"\t"$7"\t"$8"\t"$9}' plink.clumped ${inputfile}.assoc.add > ${inputfile}.assoc
                                awk 'NR==FNR{a[$1]=$0};NR>FNR{if($2 in a)print a[$2]"\t"$7"\t"$8"\t"$9}' ${inputfile}.assoc ${inputfile}.assoc.dom > ${inputfile}.dom
                                awk 'NR==FNR{a[$1]=$0};NR>FNR{if($2 in a)print a[$2]"\t"$7"\t"$8"\t"$9}'  ${inputfile}.dom ${inputfile}.assoc.rec > ${inputfile}.assoc
                                awk '{print $1}' ${inputfile}.assoc > ${inputfile}.snps
                        fi
                        echo "505"
                else
                        if [ $calculate_mod = "Origin" ];then
                                plink --bfile ${inputfile}  --pheno ${inputfile}.pheno --keep ${inputfile}.list --geno 0.05 --maf 0.05 --out ${inputfile} --linear hide-covar --covar ${sidefile} keep-pheno-on-missing-cov
                                awk '$0~/ADD/||NR==1{print $0}' ${inputfile}.assoc.linear > ${inputfile}.assoc.add
                                if [ $clumped = "False" ];then
                                        plink --bfile ${inputfile} --keep ${inputfile}.list --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump ${inputfile}.assoc.linear --clump-snp-field SNP --clump-field P
                                fi
                                awk 'NR==FNR{a[$3]=$5};NR>FNR{if(($2 in a)&&(a[$2]<0.5))print $2"\t"exp($7)"\t"$8"\t"$9}' plink.clumped ${inputfile}.assoc.add > ${inputfile}.assoc
                                awk '{print $1}' ${inputfile}.assoc > ${inputfile}.snps
                                echo "516"
                        else
                                plink --bfile ${inputfile}  --pheno ${inputfile}.pheno --keep ${inputfile}.list --geno 0.05 --maf 0.05 --out ${inputfile} --linear dominant hide-covar --covar ${sidefile} keep-pheno-on-missing-cov
                                awk '$0~/DOM/||NR==1{print $0}' ${inputfile}.assoc.linear > ${inputfile}.assoc.dom
                                plink --bfile ${inputfile}  --pheno ${inputfile}.pheno --keep ${inputfile}.list --geno 0.05 --maf 0.05 --out ${inputfile} --linear recessive hide-covar --covar ${sidefile} keep-pheno-on-missing-cov
                                awk '$0~/REC/||NR==1{print $0}' ${inputfile}.assoc.linear > ${inputfile}.assoc.rec
                                plink --bfile ${inputfile}  --pheno ${inputfile}.pheno --keep ${inputfile}.list --geno 0.05 --maf 0.05 --out ${inputfile} --linear hide-covar --covar ${sidefile} keep-pheno-on-missing-cov
                                awk '$0~/ADD/||NR==1{print $0}' ${inputfile}.assoc.linear > ${inputfile}.assoc.add
                                cat ${inputfile}.assoc.add ${inputfile}.assoc.dom ${inputfile}.assoc.rec > ${inputfile}.assoc.linear
                                if [ $clumped = "False" ];then
                                        plink --bfile ${inputfile} --keep ${inputfile}.list --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump ${inputfile}.assoc.linear --clump-snp-field SNP --clump-field P
                                        echo "527"
                                fi
                                awk 'NR==FNR{a[$3]=$5};NR>FNR{if(($2 in a)&&(a[$2]<0.5))print $2"\t"exp($7)"\t"$8"\t"$9}' plink.clumped ${inputfile}.assoc.add > ${inputfile}.assoc
                                awk 'NR==FNR{a[$1]=$0};NR>FNR{if($2 in a)print a[$2]"\t"exp($7)"\t"$8"\t"$9}' ${inputfile}.assoc ${inputfile}.assoc.dom > ${inputfile}.dom
                                awk 'NR==FNR{a[$1]=$0};NR>FNR{if($2 in a)print a[$2]"\t"exp($7)"\t"$8"\t"$9}'  ${inputfile}.dom ${inputfile}.assoc.rec > ${inputfile}.assoc
                                awk '{print $1}' ${inputfile}.assoc > ${inputfile}.snps
                        fi
                fi
                declare -A snp_site
                for snp_num_site in $(seq $minsite $takestep 0.05)
                do
                        echo "538"
                        snp_num=$(echo $input | awk -v a="$snp_num_site" 'NR>1{if($5>a){print NR-1;exit}}' plink.clumped)
                        snp_site[${#snp_site[@]}]=$snp_num
                        snp_num=$(echo $input | awk -v a="$snp_num_site" 'NR>1{if($5>a){print int((NR-2)/1000)*1000;exit}}' plink.clumped)
                        snp_site_down[${#snp_site_down[@]}]=$snp_num
                done
                echo ${snp_site[*]}
                echo ${snp_site_down[*]}
                input_num=$(echo $input | awk 'END{print NR}' ${inputfile}.snps)
                input_length=$((input_num/1000))
                split -l 1000 -d -a ${#input_length} ${inputfile}.snps ${inputfile}_
                echo "549"
                split -l 1000 -d -a ${#input_length} ${inputfile}.assoc ${inputfile}.assoc_
                echo $input_length
                echo ${#input_length}
                input_len=$((input_length-1))
        fi
        for split in $(seq 0 $input_len)
        do
        {
                split_num=`echo $split | awk '{printf("%0"'${#input_length}'"d",$0)}'`;
                echo $split_num
                echo "560"
                >${inputfile}_${split_num}.fea_few
                wc_split=$(echo $input | awk 'END{print NR}' ${inputfile}_${split_num})
                wc_ped=$(echo $input | awk 'END{print NR}' ${inputfile}.fam)
                while true
                do
                        plink --bfile ${inputfile} --keep ${inputfile}.list --extract ${inputfile}_${split_num} --recode --out ${inputfile}_${split_num}
                        plink --bfile ${inputfile} --keep ${inputfile}.list --extract ${inputfile}_${split_num} --make-bed --out ${inputfile}_${split_num}
                        wc_split_time=$(echo $input | awk 'END{print NR}' ${inputfile}_${split_num}.bim)
                        wc_ped_time=$(echo $input | awk 'END{print NR}' ${inputfile}_${split_num}.ped)
                        echo $wc_split
                        echo "571"
                        echo $wc_ped
                        break
                        if [ $wc_split_time -eq $wc_split ];then
                                if [ $wc_ped_time -eq $wc_ped ];then
                                        echo "lafhogipc;zm.mnn"
                                        break
                                fi
                        fi
                        #exit
                done
                echo "582"
                for num in $(seq 1 1000)
                do
                        echo $num
                        echo $split_num"\t"$num >> $time_point.log
                        A=$(echo $input | awk 'NR=="'$num'"{print $5}' ${inputfile}_${split_num}.bim)
                        a=$(echo $input | awk 'NR=="'$num'"{print $6}' ${inputfile}_${split_num}.bim)
                        echo $A"\t"$a"\t"$num"\t"$split_num
                        echo "593"
                        awk 'BEGIN{AA_sum=0;Aa_sum=0;aa_sum=0;dom_sum=0;rec_sum=0;het_sum=0;AA_num=0;Aa_num=0;aa_num=0;dom_num=0;rec_num=0;het_num=0;AA_std=0;Aa_std=0;aa_std=0;dom_std=0;rec_std=0;het_std=0}NR==FNR&&$0!~/-9/{if($(2*"'$num'"+5)=="'$A'"&&$(2*"'$num'"+6)=="'$A'"){AA_sum+=$6;dom_sum+=$6;het_sum+=$6;AA_num+=1;dom_num+=1;het_num+=1}else if($(2*"'$num'"+5)!=$(2*"'$num'"+6)){Aa_sum+=$6;dom_sum+=$6;rec_sum+=$6;Aa_num+=1;dom_num+=1;rec_num+=1}else if($(2*"'$num'"+5)=="'$a'"&&$(2*"'$num'"+6)=="'$a'"){aa_sum+=$6;het_sum+=$6;rec_sum+=$6;aa_num+=1;het_num+=1;rec_num+=1}}NR>FNR&&$0!~/-9/{if($(2*"'$num'"+5)=="'$A'"&&$(2*"'$num'"+6)=="'$A'"){AA_std+=($6-AA_sum/AA_num)^2;dom_std+=($6-dom_sum/dom_num)^2;het_std+=($6-het_sum/het_num)^2}else if($(2*"'$num'"+5)!=$(2*"'$num'"+6)){Aa_std+=($6-Aa_sum/Aa_num)^2;dom_std+=($6-dom_sum/dom_num)^2;rec_std+=($6-rec_sum/rec_num)^2}else if($(2*"'$num'"+5)=="'$a'"&&$(2*"'$num'"+6)=="'$a'"){aa_std+=($6-aa_sum/aa_num)^2;het_std+=($6-het_sum/het_num)^2;rec_std+=($6-rec_sum/rec_num)^2}}END{print AA_sum/AA_num"\t"AA_num"\t"sqrt(AA_std/AA_num)"\t"Aa_sum/Aa_num"\t"Aa_num"\t"sqrt(Aa_std/Aa_num)"\t"aa_sum/aa_num"\t"aa_num"\t"sqrt(aa_std/aa_num)"\t"dom_sum/dom_num"\t"dom_num"\t"sqrt(dom_std/dom_num)"\t"rec_sum/rec_num"\t"rec_num"\t"sqrt(rec_std/rec_num)"\t"het_sum/het_num"\t"het_num"\t"sqrt(het_std/het_num)"\t""'$num'"}' ${inputfile}_${split_num}.ped ${inputfile}_${split_num}.ped >> ${inputfile}_${split_num}.fea_few
                        if [ $((split*1000+num)) -eq $input_num ];then
                                break
                        fi
                done
                #rm ${inputfile}_${split_num}.ped ${inputfile}_${split_num}.log ${inputfile}_${split_num}.map
                paste ${inputfile}_${split_num}.fea_few ${inputfile}.assoc_${split_num} > ${inputfile}_${split_num}.fea
                echo "601"
                if [ $calculate_mod = "False" ];then
                        echo 101
                        #awk '{label=0;aa_ave=$7;Aa_ave=$4;AA_ave=$1;if(aa_ave>=Aa_ave&&Aa_ave>=AA_ave){if((aa_ave-Aa_ave)*3<(Aa_ave-AA_ave)){label+=4}else if((aa_ave-Aa_ave)*3>=(Aa_ave-AA_ave)&&(aa_ave-Aa_ave)<=(Aa_ave-AA_ave)*3){label+=1}else{label+=2}}else if(aa_ave>=AA_ave&&AA_ave>Aa_ave){if((aa_ave-AA_ave)>=(AA_ave-Aa_ave)){label+=2}else{label+=8}}else if(Aa_ave>aa_ave&&aa_ave>=AA_ave){if((Aa_ave-aa_ave)>(aa_ave-AA_ave)){label+=8}else{label+=4}}else if(Aa_ave>AA_ave&&AA_ave>aa_ave){if((Aa_ave-AA_ave)>(AA_ave-aa_ave)){label+=8}else{label+=2}}else if(AA_ave>aa_ave&&aa_ave>Aa_ave){if((AA_ave-aa_ave)>=(aa_ave-Aa_ave)){label+=4}else{label+=8}}else{if((AA_ave-Aa_ave)*3<(Aa_ave-aa_ave)){label+=2}else if((AA_ave-Aa_ave)*3>=(Aa_ave-aa_ave)&&(AA_ave-Aa_ave)<=(Aa_ave-aa_ave)*3){label+=1}else{label+=4}};if(label%2==1){print "add\t"log($21)"\t"sqrt($22^2)}else if(label%4==2){print "dom\t"log($23)"\t"sqrt($24^2)}else if(label%8==4){print "rec\t"(-1)*log($25)"\t"sqrt($26^2)}else{print "het\t"log(($16*$17-$17)*($5-($4*$5-$5))/(($17-($16*$17-$17))*($4*$5-$5)))"\t"sqrt((($16-$4)/sqrt($18^2/$17+$6^2/$5))^2)}}' ${inputfile}_${split_num}.fea > ${inputfile}_${split_num}.summary
                        awk '{label=0;aa_ave=$7;Aa_ave=$4;AA_ave=$1;if(aa_ave==Aa_ave){point=200}else if(aa_ave==AA_ave){point=200}else if(Aa_ave==AA_ave){point=200}else if(aa_ave>=Aa_ave&&Aa_ave>=AA_ave){if((aa_ave-Aa_ave)/(aa_ave-AA_ave)>0.5){point=((aa_ave-Aa_ave)/(aa_ave-AA_ave))*200-1000}else{point=(-1)*((aa_ave-Aa_ave)/(aa_ave-AA_ave))*200+1000}}else if(aa_ave>=AA_ave&&AA_ave>Aa_ave){point=(AA_ave-Aa_ave)/(aa_ave-AA_ave)*1000+1000}else if(Aa_ave>aa_ave&&aa_ave>=AA_ave){point=(Aa_ave-aa_ave)/(aa_ave-AA_ave)*1000+1000}else if(Aa_ave>AA_ave&&AA_ave>aa_ave){point=(Aa_ave-AA_ave)/(AA_ave-aa_ave)*1000+1000}else if(AA_ave>aa_ave&&aa_ave>Aa_ave){point=(aa_ave-Aa_ave)/(AA_ave-aa_ave)*1000+1000}else{if((AA_ave-Aa_ave)/(Aa_ave-aa_ave)>0.5){point=((AA_ave-Aa_ave)/(AA_ave-aa_ave))*200-1000}else{point=(-1)*(AA_ave-Aa_ave)/(AA_ave-aa_ave)*200+1000}};if((aa_ave-Aa_ave)^2>(Aa_ave-AA_ave)^2){label="dom"}else{label="rec"};if(label=="dom"){print "dom\t"log($24)"\t"sqrt($25^2)"\t"$26"\t"point"\t"log($21)"\t"sqrt($22^2)"\t"$23}else{print "rec\t"(-1)*log($27)"\t"sqrt($28^2)"\t"$29"\t"point"\t"log($21)"\t"sqrt($22^2)"\t"$23}}' ${inputfile}_${split_num}.fea > ${inputfile}_${split_num}.summary
                elif [ $calculate_mod = "Origin" ];then
                        echo 104
                        awk '{print "add\t"log($21)"\t"sqrt($22^2)"\t"$23}' ${inputfile}_${split_num}.fea > ${inputfile}_${split_num}.summary
                elif [ $calculate_mod = "Compare" ];then
                        echo 140
                        awk '{add_t=sqrt($22^2);dom_t=sqrt($25^2);rec_t=sqrt($28^2);if(add_t==0){point=200;print "dom\t"log($24)"\t"sqrt($25^2)"\t"$26"\t"point"\t"log($21)"\t"sqrt($22^2)"\t"$23}else if(dom_t>rec_t){point=dom_t/add_t*1000;print "dom\t"log($24)"\t"sqrt($25^2)"\t"$26"\t"point"\t"log($21)"\t"sqrt($22^2)"\t"$23}else{point=rec_t/add_t*1000;print "rec\t"(-1)*log($27)"\t"sqrt($28^2)"\t"$29"\t"point"\t"log($21)"\t"sqrt($22^2)"\t"$23}}' ${inputfile}_${split_num}.fea > ${inputfile}_${split_num}.summary
                        echo "612"
                else
                        echo 107
                        awk '{add_t=sqrt($22^2);dom_t=sqrt($25^2);rec_t=sqrt($28^2);if(log($21)==0){point=200;print "dom\t"log($24)"\t"sqrt($25^2)"\t"$26"\t"point"\t"log($21)"\t"sqrt($22^2)"\5"$23}else if(sqrt(log($24)^2)>sqrt(log($27)^2)){point=log($24)/log($21)*1000;print "dom\t"log($24)"\t"sqrt($25^2)"\t"$26"\t"point"\t"log($21)"\t"sqrt($22^2)"\5"$23}else{point=log($27)/log($21)*1000;print "rec\t"(-1)*log($27)"\t"sqrt($28^2)"\t"$29"\t"point"\t"log($21)"\t"sqrt($22^2)"\t"$23}}' ${inputfile}_${split_num}.fea > ${inputfile}_${split_num}.summary
                fi
                echo "summary OK"
                echo "split "${split_num}" summary OK" >> $time_point.log
        } &
        done
        wait
        >${inputfile}.summary_and_bim
        for split in $(seq 0 $((input_length-1)))
        do
                split_num=`echo $split | awk '{printf("%0"'${#input_length}'"d",$0)}'`;
                echo "633"
                paste ${inputfile}_${split_num}.summary ${inputfile}_${split_num}.bim >> ${inputfile}.summary_and_bim
        done
fi
if [[ $testmod == "test" ]] && [[ $validationmod == "validation" ]];then
        fam_num=$(echo $input | awk 'END{print NR}' $inputfile.fam)
        threshold_t=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c "import t_to_p;print(t_to_p.t_to_p("${threshold_t},${fam_num}"))"`
        echo "threshold_t\t"$threshold_t
elif [[ $testmod != "test" ]] && [[ $validationmod == "validation" ]];then
        fam_num=$(echo $input | awk 'END{print NR*9/10}' $inputfile.fam)
        threshold_t=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c "import t_to_p;print(t_to_p.t_to_p("${threshold_t},${fam_num}"))"`
        echo "threshold_t\t"$threshold_t
        echo "296"
elif [[ $testmod == "test" ]] && [[ $validationmod != "validation" ]];then
        fam_num=$(echo $input | awk 'END{print NR*9/10}' $inputfile.fam)
        threshold_t=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c "import t_to_p;print(t_to_p.t_to_p("${threshold_t},${fam_num}"))"`
        echo "threshold_t\t"$threshold_t
else
        fam_num=$(echo $input | awk 'END{print NR*4/5}' $inputfile.fam)
        fam_num_validation=$(echo $input | awk 'END{print NR*9/10}' $inputfile.fam)
        threshold_t=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c "import t_to_p;print(t_to_p.t_to_p("${threshold_t},${fam_num}"))"`
        echo "threshold_t\t"$threshold_t
fi
if [[ $testmod != "test" ]] && [[ $validationmod == "validation" ]];then
        awk -v line_a=$fam_num 'NR>line_a&&NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${inputfile}.fam ${keep_file} >> ${testfile}.list
        awk -v line_a=$fam_num 'NR>1&&NR<=line_a{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${inputfile}.fam ${keep_file} >> ${inputfile}.list
        awk 'NR>1&&NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${validationfile}.fam ${keep_file} >> ${validationfile}.list
elif [[ $testmod == "test" ]] && [[ $validationmod != "validation" ]];then
        awk -v line_a=$fam_num 'NR>line_a&&NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${inputfile}.fam ${keep_file} >> ${validationfile}.list
        awk -v line_a=$fam_num 'NR>1&&NR<=line_a&&NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${inputfile}.fam ${keep_file} >> ${inputfile}.list
        awk 'NR>1&&NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${testfile}.fam ${keep_file} >> ${testfile}.list
elif [[ $testmod == "test" ]] && [[ $validationmod == "validation" ]];then
        echo "329"
        awk 'NR>1&&NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${validationfile}.fam ${keep_file} >> ${validationfile}.list
        awk 'NR>1&&NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${inputfile}.fam ${keep_file} >> ${inputfile}.list
        awk 'NR>1&&NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${testfile}.fam ${keep_file} >> ${testfile}.list
else
        awk -v line_a=$fam_num_validation 'NR>line_a&&NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${inputfile}.fam ${keep_file} >> ${testfile}.list
        awk -v line_a=$fam_num 'NR>1&&NR<=line_a&&NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${inputfile}.fam ${keep_file} >> ${inputfile}.list
        awk -v line_a=$fam_num -v line_b=$fam_num_validation 'NR>line_a&&NR<=line_b&&NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $1"\t"$2}}' ${inputfile}.fam ${keep_file} >> ${validationfile}.list
fi
comm -12 <(sort ${inputfile}.list|uniq) <(sort ${inputfile}.list|uniq) > time
mv time ${inputfile}.list
comm -12 <(sort ${validationfile}.list|uniq) <(sort ${validationfile}.list|uniq) > time
mv time ${validationfile}.list
comm -12 <(sort ${testfile}.list|uniq) <(sort ${testfile}.list|uniq) > time
mv time ${testfile}.list
if [ $calculate_mod = "Origin" ];then
        sort -g -k 4 ${inputfile}.summary_and_bim -o time
else
        sort -g -k 8 ${inputfile}.summary_and_bim -o time
fi
mv time ${inputfile}.summary_and_bim
if [ $calculate_mod = "Origin" ];then
        awk '{print $6}' ${inputfile}.summary_and_bim > ${inputfile}.snps
else
        awk '{print $10}' ${inputfile}.summary_and_bim > ${inputfile}.snps
fi
if [ $group = "True" ] ;then
        for inputfile_order in $(seq 1 $inputfile_num)
        do
                plink --bfile ${inputfile_order} --keep ${inputfile}.list --extract ${inputfile}.snps --recode --out ${inputfile_order}_${split_num}
                echo "644"
                plink --bfile ${inputfile_order} --keep ${inputfile}.list --extract ${inputfile}.snps --make-bed --out ${inputfile_order}_${split_num}
                echo ${inputfile_order}_${split_num} > time
        done
        head -1 time > time2
        sed -i '1d' time
        plink --keep ${inputfile}.list --noweb --bfile time2 --merge-list time --make-bed --out ${inputfile}
        plink --keep ${inputfile}.list --noweb --bfile time2 --merge-list time --recode --out ${inputfile}
fi
input_num=$(echo $input | awk 'END{print NR}' ${inputfile}.summary_and_bim)
input_length=$((input_num/1000))
echo "655"
split -l 1000 -d -a ${#input_length} ${inputfile}.summary_and_bim ${inputfile}.summary_and_bim_
for threshold in $(seq $begin_site $step $end_site)
do
        awk '{print "0"}' ${testfile}.list > time_${threshold}.best
        awk '{print "0"}' ${testfile}.list > time.best
        echo "time_"${threshold}".best"
        head -2 time_${threshold}.best
        echo "time.best"
        head -2 time.best
        awk '{print "0"}' ${validationfile}.list > time_${threshold}_validation.best
        awk '{print "0"}' ${validationfile}.list > time_validation.best
        echo "time_"${threshold}"_validation.best"
        head -2 time_${threshold}_validation.best
        echo "time_validation.best"
        head -2 time_validation.best
done
declare -A old_relation
echo "666"
for threshold in $(seq $begin_site $step $end_site)
do
        >$threshold.log
        old_relation[$threshold]=0
done
wc -l ${inputfile}.snps
head ${inputfile}.snps
input_num=$(echo $input | awk 'END{print NR}' ${inputfile}.snps) input_length=$((input_num/1000))
split -l 1000 -d -a ${#input_length} ${inputfile}.snps ${inputfile}_
declare -A snp_site
for snp_num_site in $(seq $minsite $takestep 0.05)
do
        echo $snp_num_site
        if [ $calculate_mod = "Origin" ];then
                snp_num=$(echo $input | awk -v a="$snp_num_site" 'NR>1{if($8>a){print NR-1;exit}}' ${inputfile}.summary_and_bim)
                echo $snp_num
                snp_site[${#snp_site[@]}]=$snp_num
                snp_num=$(echo $input | awk -v a="$snp_num_site" 'NR>1{if($8>a){print int((NR-2)/1000)*1000;exit}}' ${inputfile}.summary_and_bim)
                echo "680"
                echo $snp_num
                snp_site_down[${#snp_site_down[@]}]=$snp_num
        else
                snp_num=$(echo $input | awk -v a="$snp_num_site" 'NR>1{if($4>a){print NR-1;exit}}' ${inputfile}.summary_and_bim)
                echo $snp_num
                snp_site[${#snp_site[@]}]=$snp_num
                snp_num=$(echo $input | awk -v a="$snp_num_site" 'NR>1{if($4>a){print int((NR-2)/1000)*1000;exit}}' ${inputfile}.summary_and_bim)
                echo "680"
                echo $snp_num
                snp_site_down[${#snp_site_down[@]}]=$snp_num
        fi
done
echo ${snp_site[*]}
echo ${snp_site_down[*]}
echo $input_length
echo ${#input_length}
input_len=$((input_length-1))
if [ $prsice = "False" ];then
        split -l 1000 -d -a ${#input_length} ${inputfile}.snps ${inputfile}_
        echo "688"
        echo $input_len
        for spli_num in $(seq 0 $input_len)
        do
                echo "sort OK"
                echo "sort OK" >> $time_point.log
                tail ${inputfile}_${split_num}
                split_num=`echo $spli_num | awk '{printf("%0"'${#input_length}'"d",$0)}'`;
                if [[ $testmod == "test" ]] && [[ $validationmod == "validation" ]];then
                        plink --bfile ${testfile} --keep ${testfile}.list --extract ${inputfile}_${split_num} --recode --out ${testfile}_${split_num}
                        plink --bfile ${testfile} --keep ${testfile}.list --extract ${inputfile}_${split_num} --make-bed --out ${testfile}_${split_num}
                        plink --bfile ${validationfile} --keep ${validationfile}.list --extract ${inputfile}_${split_num} --recode --out ${validationfile}_${split_num}
                        plink --bfile ${validationfile} --keep ${validationfile}.list --extract ${inputfile}_${split_num} --make-bed --out ${validationfile}_${split_num}
                elif [[ $testmod != "test" ]] && [[ $validationmod == "validation" ]];then
                        plink --bfile ${inputfile} --keep ${testfile}.list --extract ${inputfile}_${split_num} --recode --out ${testfile}_${split_num}
                        plink --bfile ${inputfile} --keep ${testfile}.list --extract ${inputfile}_${split_num} --make-bed --out ${testfile}_${split_num}
                        plink --bfile ${validationfile} --keep ${validationfile}.list --extract ${inputfile}_${split_num} --recode --out ${validationfile}_${split_num}
                        plink --bfile ${validationfile} --keep ${validationfile}.list --extract ${inputfile}_${split_num} --make-bed --out ${validationfile}_${split_num}
                elif [[ $testmod == "test" ]] && [[ $validationmod != "validation" ]];then
                        plink --bfile ${testfile} --keep ${testfile}.list --extract ${inputfile}_${split_num} --recode --out ${testfile}_${split_num}
                        plink --bfile ${testfile} --keep ${testfile}.list --extract ${inputfile}_${split_num} --make-bed --out ${testfile}_${split_num}
                        plink --bfile ${inputfile} --keep ${validationfile}.list --extract ${inputfile}_${split_num} --recode --out ${validationfile}_${split_num}
                        plink --bfile ${inputfile} --keep ${validationfile}.list --extract ${inputfile}_${split_num} --make-bed --out ${validationfile}_${split_num}
                else
                        head ${testfile}.list
                        echo "${testfile}.list"
                        head ${inputfile}.fam
                        echo "${inputfile}.fam"
                        plink --bfile ${inputfile} --keep ${testfile}.list --extract ${inputfile}_${split_num} --recode --out ${testfile}_${split_num}
                        plink --bfile ${inputfile} --keep ${testfile}.list --extract ${inputfile}_${split_num} --make-bed --out ${testfile}_${split_num}
                        plink --bfile ${inputfile} --keep ${validationfile}.list --extract ${inputfile}_${split_num} --recode --out ${validationfile}_${split_num}
                        plink --bfile ${inputfile} --keep ${validationfile}.list --extract ${inputfile}_${split_num} --make-bed --out ${validationfile}_${split_num}
                fi

                echo "454"
                echo "${snp_site_down[@]}" | grep -wq $((1000*10#$split_num)) && symbol=1 || echo "457"

                if [ $symbol = 1 ];then
                        for inputnum in $(seq 1 1000)
                        do
                                if [ $calculate_mod = "Origin" ];then
                                        num=$(echo $input | awk 'NR=="'$inputnum'"{a=$6}NR>FNR{if($2==a){print FNR}}' ${inputfile}.summary_and_bim_${split_num} ${testfile}_${split_num}.map)
                                        echo "730"
                                        awk 'NR=="'$inputnum'"{c=$9;beta=$2;type=$1;t=$3}ARGIND==2{human=FNR;if($(2*"'$num'"+5)==c&&$(2*"'$num'"+6)==c&&$(2*"'$num'"+5)!=0){sum[FNR]+=beta*2}else if($(2*"'$num'"+5)!=$(2*"'$num'"+6)||$(2*"'$num'"+5)=="0"){sum[FNR]+=beta}else{sum[FNR]+=0}}ARGIND==3{print sum[FNR]+$1}' ${inputfile}.summary_and_bim_${split_num} ${testfile}_${split_num}.ped time.best > new_time.best
                                        mv new_time.best time.best
                                        awk 'NR=="'$inputnum'"{c=$9;beta=$2;type=$1;t=$3}NR>FNR{human=FNR;if($(2*"'$num'"+5)==c&&$(2*"'$num'"+6)==c&&$(2*"'$num'"+5)!=0){sum[FNR]+=beta*2}else if($(2*"'$num'"+5)!=$(2*"'$num'"+6)||$(2*"'$num'"+5)=="0"){sum[FNR]+=beta}else{sum[FNR]+=0}}ARGIND==3{print sum[FNR]+$1}' ${inputfile}.summary_and_bim_${split_num} ${validationfile}_${split_num}.ped time_validation.best > new_time_validation.best
                                        mv new_time_validation.best time_validation.best
                                elif [ $merge = "False" ];then
                                        num=$(echo $input | awk 'NR=="'$inputnum'"{a=$10}NR>FNR{if($2==a){print FNR}}' ${inputfile}.summary_and_bim_${split_num} ${testfile}_${split_num}.map)
                                        echo "730"
                                        for threshold in $(seq $begin_site $step $end_site)
                                        do
                                                {
                                                        awk -v a="$threshold" 'NR=="'$inputnum'"{c=$13;beta_type=$2;type=$1;beta_add=$6;point=$5}ARGIND==2{human=FNR;if(point<a){if($(2*"'$num'"+5)==c&&$(2*"'$num'"+6)==c&&$(2*"'$num'"+5)!=0){sum[FNR]+=beta_add*2}else if($(2*"'$num'"+5)!=$(2*"'$num'"+6)||$(2*"'$num'"+5)=="0"){sum[FNR]+=beta_add}else{sum[FNR]+=0}}else if(type=="dom"){if(($(2*"'$num'"+5)==c&&$(2*"'$num'"+6)==c&&$(2*"'$num'"+5)!=0)||($(2*"'$num'"+5)!=$(2*"'$num'"+6))){sum[FNR]+=beta_type*2}else if($(2*"'$num'"+5)=="0"){sum[FNR]+=beta_type}else{sum[FNR]+=0}}else if(type=="rec"){if(($(2*"'$num'"+5)==c&&$(2*"'$num'"+6)==c&&$(2*"'$num'"+5)!=0)){sum[FNR]-=beta_type*2}else if($(2*"'$num'"+5)=="0"){sum[FNR]-=beta_type}else{sum[FNR]+=0}}}ARGIND==3{print sum[FNR]+$1}' ${inputfile}.summary_and_bim_${split_num} ${testfile}_${split_num}.ped time_${threshold}.best > new_time_${threshold}.best
                                                        echo "741"
                                                        mv new_time_${threshold}.best time_${threshold}.best
                                                        awk -v a="$threshold" 'NR=="'$inputnum'"{c=$13;beta_type=$2;type=$1;beta_add=$6;point=$5}ARGIND==2{human=FNR;if(point<a){if($(2*"'$num'"+5)==c&&$(2*"'$num'"+6)==c&&$(2*"'$num'"+5)!=0){sum[FNR]+=beta_add*2}else if($(2*"'$num'"+5)!=$(2*"'$num'"+6)||$(2*"'$num'"+5)=="0"){sum[FNR]+=beta_add}else{sum[FNR]+=0}}else if(type=="dom"){if(($(2*"'$num'"+5)==c&&$(2*"'$num'"+6)==c&&$(2*"'$num'"+5)!=0)||($(2*"'$num'"+5)!=$(2*"'$num'"+6))){sum[FNR]+=beta_type*2}else if($(2*"'$num'"+5)=="0"){sum[FNR]+=beta_type}else{sum[FNR]+=0}}else if(type=="rec"){if(($(2*"'$num'"+5)==c&&$(2*"'$num'"+6)==c&&$(2*"'$num'"+5)!=0)){sum[FNR]-=beta_type*2}else if($(2*"'$num'"+5)=="0"){sum[FNR]-=beta_type}else{sum[FNR]+=0}}}ARGIND==3{print sum[FNR]+$1}' ${inputfile}.summary_and_bim_${split_num} ${validationfile}_${split_num}.ped time_${threshold}_validation.best > new_time_${threshold}_validation.best
                                                        mv new_time_${threshold}_validation.best time_${threshold}_validation.best
                                                } &
                                        done
                                        wait
                                        if [ $with_origin = "True" ];then
                                                awk 'NR=="'$inputnum'"{c=$13;beta=$2;type=$1;t=$3}ARGIND==2{human=FNR;if($(2*"'$num'"+5)==c&&$(2*"'$num'"+6)==c&&$(2*"'$num'"+5)!=0){sum[FNR]+=beta*2}else if($(2*"'$num'"+5)!=$(2*"'$num'"+6)||$(2*"'$num'"+5)=="0"){sum[FNR]+=beta}else{sum[FNR]+=0}}ARGIND==3{print sum[FNR]+$1}' ${inputfile}.summary_and_bim_${split_num} ${testfile}_${split_num}.ped time.best > new_time.best
                                                mv new_time.best time.best
                                                awk 'NR=="'$inputnum'"{c=$13;beta=$2;type=$1;t=$3}NR>FNR{human=FNR;if($(2*"'$num'"+5)==c&&$(2*"'$num'"+6)==c&&$(2*"'$num'"+5)!=0){sum[FNR]+=beta*2}else if($(2*"'$num'"+5)!=$(2*"'$num'"+6)||$(2*"'$num'"+5)=="0"){sum[FNR]+=beta}else{sum[FNR]+=0}}ARGIND==3{print sum[FNR]+$1}' ${inputfile}.summary_and_bim_${split_num} ${validationfile}_${split_num}.ped time_validation.best > new_time_validation.best
                                                mv new_time_validation.best time_validation.best
                                        fi
                                else
                                        num=$(echo $input | awk 'NR=="'$inputnum'"{a=$10}NR>FNR{if($2==a){print FNR}}' ${inputfile}.summary_and_bim_${split_num} ${testfile}_${split_num}.map)
                                        echo "730"
                                        for threshold in $(seq $begin_site $step $end_site)
                                        do
                                                {
                                                        echo "752"
                                                        awk -v a="$threshold" 'function max(a,b){if(a>b){return a}else{return b}};function min(a,b){if(a<b){return a}else{return b}};NR=="'$inputnum'"{c=$13;beta_type=$2;type=$1;beta_add=$6;point=$5}ARGIND==2{human=FNR;if($(2*"'$num'"+5)==c&&$(2*"'$num'"+6)==c&&$(2*"'$num'"+5)!=0){sum[FNR]+=beta_add*2*max(0,a-point)}else if($(2*"'$num'"+5)!=$(2*"'$num'"+6)||$(2*"'$num'"+5)=="0"){sum[FNR]+=beta_add*max(0,a-point)}else{sum[FNR]+=0}if(type=="dom"){if(($(2*"'$num'"+5)==c&&$(2*"'$num'"+6)==c&&$(2*"'$num'"+5)!=0)||($(2*"'$num'"+5)!=$(2*"'$num'"+6))){sum[FNR]+=beta_type*2*min(a,point)}else if($(2*"'$num'"+5)=="0"){sum[FNR]+=beta_type*min(a,point)}else{sum[FNR]+=0}}else if(type=="rec"){if(($(2*"'$num'"+5)==c&&$(2*"'$num'"+6)==c&&$(2*"'$num'"+5)!=0)){sum[FNR]-=beta_type*2*min(a,point)}else if($(2*"'$num'"+5)=="0"){sum[FNR]-=beta_type*min(a,point)}else{sum[FNR]+=0}}}ARGIND==3{print sum[FNR]+$1}' ${inputfile}.summary_and_bim_${split_num} ${testfile}_${split_num}.ped time_${threshold}.best > new_time_${threshold}.best
                                                        echo ${inputfile}".summary_and_bim_"${split_num}
                                                        head -2 ${inputfile}.summary_and_bim_${split_num}
                                                        echo ${testfile}"_"${split_num}".ped"
                                                        head -2 ${testfile}_${split_num}.ped
                                                        echo "time_"${threshold}".best"
                                                        head -2 time_${threshold}.best
                                                        mv new_time_${threshold}.best time_${threshold}.best
                                                        awk -v a="$threshold" 'function max(a,b){if(a>b){return a}else{return b}};function min(a,b){if(a<b){return a}else{return b}};NR=="'$inputnum'"{c=$13;beta_type=$2;type=$1;beta_add=$6;point=$5}ARGIND==2{human=FNR;if($(2*"'$num'"+5)==c&&$(2*"'$num'"+6)==c&&$(2*"'$num'"+5)!=0){sum[FNR]+=beta_add*2*max(0,a-point)}else if($(2*"'$num'"+5)!=$(2*"'$num'"+6)||$(2*"'$num'"+5)=="0"){sum[FNR]+=beta_add*max(0,a-point)}else{sum[FNR]+=0}if(type=="dom"){if(($(2*"'$num'"+5)==c&&$(2*"'$num'"+6)==c&&$(2*"'$num'"+5)!=0)||($(2*"'$num'"+5)!=$(2*"'$num'"+6))){sum[FNR]+=beta_type*2*min(a,point)}else if($(2*"'$num'"+5)=="0"){sum[FNR]+=beta_type*min(a,point)}else{sum[FNR]+=0}}else if(type=="rec"){if(($(2*"'$num'"+5)==c&&$(2*"'$num'"+6)==c&&$(2*"'$num'"+5)!=0)){sum[FNR]-=beta_type*2*min(a,point)}else if($(2*"'$num'"+5)=="0"){sum[FNR]-=beta_type*min(a,point)}else{sum[FNR]+=0}}}ARGIND==3{print sum[FNR]+$1}' ${inputfile}.summary_and_bim_${split_num} ${validationfile}_${split_num}.ped time_${threshold}_validation.best > new_time_${threshold}_validation.best
                                                        echo ${validationfile}"_"${split_num}".ped"
                                                        head -2 ${validationfile}_${split_num}.ped
                                                        echo "time_"${threshold}"_validation.best"
                                                        head -2 time_${threshold}_validation.best
                                                        mv new_time_${threshold}_validation.best time_${threshold}_validation.best
                                                } &
                                        done
                                        wait
                                        if [ $with_origin = "True" ];then
                                                awk 'NR=="'$inputnum'"{c=$13;beta=$2;type=$1;t=$3}ARGIND==2{human=FNR;if($(2*"'$num'"+5)==c&&$(2*"'$num'"+6)==c&&$(2*"'$num'"+5)!=0){sum[FNR]+=beta*2}else if($(2*"'$num'"+5)!=$(2*"'$num'"+6)||$(2*"'$num'"+5)=="0"){sum[FNR]+=beta}else{sum[FNR]+=0}}ARGIND==3{print sum[FNR]+$1}' ${inputfile}.summary_and_bim_${split_num} ${testfile}_${split_num}.ped time.best > new_time.best
                                                mv new_time.best time.best
                                                awk 'NR=="'$inputnum'"{c=$13;beta=$2;type=$1;t=$3}NR>FNR{human=FNR;if($(2*"'$num'"+5)==c&&$(2*"'$num'"+6)==c&&$(2*"'$num'"+5)!=0){sum[FNR]+=beta*2}else if($(2*"'$num'"+5)!=$(2*"'$num'"+6)||$(2*"'$num'"+5)=="0"){sum[FNR]+=beta}else{sum[FNR]+=0}}ARGIND==3{print sum[FNR]+$1}' ${inputfile}.summary_and_bim_${split_num} ${validationfile}_${split_num}.ped time_validation.best > new_time_validation.best
                                                mv new_time_validation.best time_validation.best
                                        fi
                                fi

                                echo $((inputnum+1000*10#$split_num))"\t749"
                                echo $((inputnum))"\t750"
                                echo "${snp_site[@]}" | grep -wq $((inputnum+1000*10#$split_num)) && echo "495" || continue
                                echo "763"

                                if [ $calculate_mod != "Origin" ];then
                                        for threshold in $(seq $begin_site $step $end_site)
                                        do
                                        {
                                                #awk 'NR==FNR{a[NR]=$1}NR>FNR{print a[FNR]+$1}' ${num}_${split_num}_${threshold}.best all_${threshold}.best > time_${threshold}.best
                                                #awk 'NR==FNR{a[NR]=$1}NR>FNR{print a[FNR]+$1}' ${num}_${split_num}_${threshold}_validation.best all_${threshold}_validation.best > time_${threshold}_validation.best
                                                paste time_${threshold}.best ${testfile}_${split_num}.fam > time_${threshold}.spearman
                                                echo "time_'${threshold}'.best"
                                                head -2 time_${threshold}.best
                                                echo ${testfile}"_"${split_num}".fam"
                                                head -2 ${testfile}_${split_num}.fam
                                                paste time_${threshold}_validation.best ${validationfile}_${split_num}.fam > time_${threshold}_validation.spearman
                                                echo "time_'${threshold}'_validation.best"
                                                head -2 time_${threshold}_validation.best
                                                echo ${validationfile}"_"${split_num}".fam"
                                                head -2 ${validationfile}_${split_num}.fam
                                                sed -i '/-9/d' time_${threshold}.spearman
                                                echo "774"
                                                sed -i '/-9/d' time_${threshold}_validation.spearman
                                                if [ $logistic = "False" ];then
                                                        old_relation_validation_spearman=-1000
                                                        old_relation_validation_pearson=-1000
                                                        old_relation_validation_r2=-1000
                                                        relation_spearman=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import spearman;print(spearman.ReadData("time_'${threshold}'.spearman"))'`
                                                        relation_validation_spearman=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import spearman;print(spearman.ReadData("time_'${threshold}'_validation.spearman"))'`
                                                        relation_r2=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import r2;print(r2.ReadData("time_'${threshold}'.spearman"))'`
                                                        relation_validation_r2=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import r2;print(r2.ReadData("time_'${threshold}'_validation.spearman"))'`
                                                        relation_pearson=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import pearson;print(pearson.ReadData("time_'${threshold}'.spearman"))'`
                                                        relation_validation_pearson=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import pearson;print(pearson.ReadData("time_'${threshold}'_validation.spearman"))'`
                                                        echo "the "$((inputnum+1000*10#$split_num))" snp by the "$threshold" threshold got results as spearman values: test set as "$relation_spearman " and validation set as "$relation_validation_spearman"; r2 values: test set as "$relation_r2 " and validation set as "$relation_validation_r2"; pearson values: test set as "$relation_pearson " and validation set as "$relation_validation_pearson
                                                        echo "the "$((inputnum+1000*10#$split_num))" snp by the "$threshold" threshold got results as spearman values: test set as "$relation_spearman " and validation set as "$relation_validation_spearman"; r2 values: test set as "$relation_r2 " and validation set as "$relation_validation_r2"; pearson values: test set as "$relation_pearson " and validation set as "$relation_validation_pearson >> $time_point.log
                                                        echo "785"
                                                        echo "the "$((inputnum+1000*10#$split_num))" snp by the "$threshold" threshold got results as spearman values: test set as "$relation_spearman " and validation set as "$relation_validation_spearman"; r2 values: test set as "$relation_r2 " and validation set as "$relation_validation_r2"; pearson values: test set as "$relation_pearson " and validation set as "$relation_validation_pearson >> ${threshold}.log
                                                        sort -g -k 2 ${threshold}.log -o ${threshold}.log
                                                else
                                                        old_relation_validation_bce=-1000
                                                        old_relation_validation_r2=-1000
                                                        relation_bce=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import bce;print(bce.ReadData("time_'${threshold}'.spearman"))'`
                                                        echo "time_'${threshold}'.spearman"
                                                        head -2 time_${threshold}.spearman
                                                        relation_validation_bce=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import bce;print(bce.ReadData("time_'${threshold}'_validation.spearman"))'`
                                                        echo "time_'${threshold}'_validation.spearman"
                                                        head -2 time_${threshold}_validation.spearman
                                                        relation_r2=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import r2;print(r2.ReadData("time_'${threshold}'.spearman"))'`
                                                        echo "time_'${threshold}'.spearman"
                                                        head -2 time_${threshold}.spearman
                                                        relation_validation_r2=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import r2;print(r2.ReadData("time_'${threshold}'_validation.spearman"))'`
                                                        echo "time_'${threshold}'_validation.spearman"
                                                        head -2 time_${threshold}_validation.spearman
                                                        echo "the "$((inputnum+1000*10#$split_num))" snp by the "$threshold" threshold got results as AUCROC values: test set as "$relation_bce " and validation set as "$relation_validation_bce"; r2 values: test set as "$relation_r2 " and validation set as "$relation_validation_r2
                                                        echo "the "$((inputnum+1000*10#$split_num))" snp by the "$threshold" threshold got results as AUCROC values: test set as "$relation_bce " and validation set as "$relation_validation_bce"; r2 values: test set as "$relation_r2 " and validation set as "$relation_validation_r2 >> $time_point.log
                                                        echo "the "$((inputnum+1000*10#$split_num))" snp by the "$threshold" threshold got results as AUCROC values: test set as "$relation_bce " and validation set as "$relation_validation_bce"; r2 values: test set as "$relation_r2 " and validation set as "$relation_validation_r2 >> ${threshold}.log
                                                        sort -g -k 2 $threshold.log -o $threshold.log
                                                fi
                                                echo "796"
                                                #mv time_${threshold}.best all_${threshold}.best
                                                #mv time_${threshold}_validation.best all_${threshold}_validation.best
                                        } &
                                        done
                                        wait
                                fi
                                if [[ $calculate_mod == "Origin" ]] || [[ $with_origin == "True" ]];then
                                        #awk 'NR==FNR{a[NR]=$1}NR>FNR{print a[FNR]+$1}' ${num}_${split_num}.best all.best > time.best
                                        #awk 'NR==FNR{a[NR]=$1}NR>FNR{print a[FNR]+$1}' ${num}_${split_num}_validation.best all_validation.best > time_validation.best
                                        paste time.best ${testfile}_${split_num}.fam > time.spearman
                                        paste time_validation.best ${validationfile}_${split_num}.fam > time_validation.spearman
                                        echo "807"
                                        sed -i '/-9/d' time.spearman
                                        sed -i '/-9/d' time_validation.spearman
                                        if [ $logistic = "False" ];then
                                                old_relation_validation_spearman=-1000
                                                old_relation_validation_pearson=-1000
                                                old_relation_validation_r2=-1000
                                                relation_spearman=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import spearman;print(spearman.ReadData("time.spearman"))'`
                                                relation_validation_spearman=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import spearman;print(spearman.ReadData("time_validation.spearman"))'`
                                                relation_r2=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import r2;print(r2.ReadData("time.spearman"))'`
                                                relation_validation_r2=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import r2;print(r2.ReadData("time_validation.spearman"))'`
                                                relation_pearson=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import pearson;print(pearson.ReadData("time.spearman"))'`
                                                relation_validation_pearson=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import pearson;print(pearson.ReadData("time_validation.spearman"))'`
                                                echo "the "$((inputnum+1000*10#$split_num))" snp got results as spearman values: test set as "$relation_spearman " and validation set as "$relation_validation_spearman"; r2 values: test set as "$relation_r2 " and validation set as "$relation_validation_r2"; pearson values: test set as "$relation_pearson " and validation set as "$relation_validation_pearson
                                                echo "818"
                                                echo "the "$((inputnum+1000*10#$split_num))" snp threshold got results as spearman values: test set as "$relation_spearman " and validation set as "$relation_validation_spearman"; r2 values: test set as "$relation_r2 " and validation set as "$relation_validation_r2"; pearson values: test set as "$relation_pearson " and validation set as "$relation_validation_pearson >> ${time_point}.log
                                                sort -g -k 2 ${time_point}.log -o ${time_point}.log
                                        else
                                                old_relation_validation_bce=-1000
                                                old_relation_validation_r2=-1000
                                                relation_bce=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import bce;print(bce.ReadData("time.spearman"))'`
                                                relation_validation_bce=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import bce;print(bce.ReadData("time_validation.spearman"))'`
                                                relation_r2=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import r2;print(r2.ReadData("time.spearman"))'`
                                                relation_validation_r2=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import r2;print(r2.ReadData("time_validation.spearman"))'`
                                                echo "the "$((inputnum+1000*10#$split_num))" snp got results as AUCROC values: test set as "$relation_bce " and validation set as "$relation_validation_bce
                                                echo "the "$((inputnum+1000*10#$split_num))" snp got results as AUCROC values: test set as "$relation_bce " and validation set as "$relation_validation_bce >> ${time_point}.log
                                                echo "the "$((inputnum+1000*10#$split_num))" snp got results as AUCROC values: test set as "$relation_bce " and validation set as "$relation_validation_bce"; r2 values: test set as "$relation_r2 " and validation set as "$relation_validation_r2
                                                echo "the "$((inputnum+1000*10#$split_num))" snp got results as AUCROC values: test set as "$relation_bce " and validation set as "$relation_validation_bce"; r2 values: test set as "$relation_r2 " and validation set as "$relation_validation_r2 >> $time_point.log
                                                echo "the "$((inputnum+1000*10#$split_num))" snp got results as AUCROC values: test set as "$relation_bce " and validation set as "$relation_validation_bce"; r2 values: test set as "$relation_r2 " and validation set as "$relation_validation_r2 >> ${threshold}.log
                                                sort -g -k 2 ${time_point}.log -o ${time_point}.log
                                        fi
                                fi
                                echo "540"
                        done
                        if [ $calculate_mod != "Origin" ];then
                                sort -g -k 2 ${threshold}.log -o ${threshold}.log
                                /home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import draw;draw("$begin_site $step $end_site")'
                                echo "840"
                                sort -r -g -k 21 ${threshold}.log -o ${threshold}.log
                        else
                                sort -g -k 2 ${time_point}.log -o ${time_point}.log
                                /home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import draw;draw_origin("$time_point")'
                                sort -r -g -k 21 ${time_point}.log -o ${time_point}.log
                        fi
                fi

                if [ $symbol = 0 ];then
                        echo "852"
                        if [ $calculate_mod = "Origin" ];then
                                awk 'NR==FNR{c[FNR]=$9;beta_add[FNR]=$2;type_add[FNR]=$1;t_add[FNR]=$3;rs[FNR]=$6}ARGIND==2{for(snp=1;snp<=1000;snp++){if($2==rs[snp]){snp_site[FNR]=snp}}}ARGIND==3{human=FNR;for(snp=1;snp<=1000;snp++){if($(2*snp+5)==c[snp_site[snp]]&&$(2*snp+6)==c[snp_site[snp]]&&$(2*snp+5)!=0){sum[FNR]+=beta_add[snp_site[snp]]*2}else if($(2*snp+5)!=$(2*snp+6)||$(2*snp+5)==0){sum[FNR]+=beta_add[snp_site[snp]]}else{sum[FNR]+=0}}}ARGIND==4{print sum[FNR]+$1}' ${inputfile}.summary_and_bim_${split_num} ${testfile}_${split_num}.map ${testfile}_${split_num}.ped time.best > new_time.best
                                mv new_time.best time.best
                                awk 'NR==FNR{c[FNR]=$9;beta_add[FNR]=$2;type_add[FNR]=$1;t_add[FNR]=$3;rs[FNR]=$6}ARGIND==2{for(snp=1;snp<=1000;snp++){if($2==rs[snp]){snp_site[FNR]=snp}}}ARGIND==3{human=FNR;for(snp=1;snp<=1000;snp++){if($(2*snp+5)==c[snp_site[snp]]&&$(2*snp+6)==c[snp_site[snp]]&&$(2*snp+5)!=0){sum[FNR]+=beta_add[snp_site[snp]]*2}else if($(2*snp+5)!=$(2*snp+6)||$(2*snp+5)==0){sum[FNR]+=beta_add[snp_site[snp]]}else{sum[FNR]+=0}}}ARGIND==4{print sum[FNR]+$1}' ${inputfile}.summary_and_bim_${split_num} ${validationfile}_${split_num}.map ${validationfile}_${split_num}.ped time_validation.best > new_time_validation.best
                                mv new_time_validation.best time_validation.best
                        elif [ $merge = "False" ];then
                                for threshold in $(seq $begin_site $step $end_site)
                                do
                                        {
                                                awk -v a="$threshold" 'NR==FNR{c[FNR]=$13;beta_add[FNR]=$6;type[FNR]=$1;beta_type[FNR]=$2;point[FNR]=$5;rs[FNR]=$10}ARGIND==2{for(snp=1;snp<=1000;snp++){if($2==rs[snp]){snp_site[FNR]=snp}}}ARGIND==3{human=FNR;for(snp=1;snp<=1000;snp++){if(point[snp_site[snp]]<a){if($(2*snp+5)==c[snp_site[snp]]&&$(2*snp+6)==c[snp_site[snp]]&&$(2*snp+5)!=0){sum[FNR]+=beta_add[snp_site[snp]]*2}else if($(2*snp+5)!=$(2*snp+6)||$(2*snp+5)==0){sum[FNR]+=beta_add[snp_site[snp]]}else{sum[FNR]+=0}}else if(type[snp_site[snp]]=="dom"){if(($(2*snp+5)==c[snp_site[snp]]&&$(2*snp+6)==c[snp_site[snp]]&&$(2*snp+5)!=0)||($(2*snp+5)!=$(2*snp+6))){sum[FNR]+=beta_type[snp_site[snp]]*2}else if($(2*snp+5)==0){sum[FNR]+=beta_type[snp_site[snp]]}else{sum[FNR]+=0}}else{if(($(2*snp+5)==c[snp_site[snp]]&&$(2*snp+6)==c[snp_site[snp]]&&$(2*snp+5)!=0)){sum[FNR]-=beta_type[snp_site[snp]]*2}else if($(2*snp+5)==0){sum[FNR]-=beta_type[snp_site[snp]]}else{sum[FNR]+=0}}}}ARGIND==4{print sum[FNR]+$1}' ${inputfile}.summary_and_bim_${split_num} ${testfile}_${split_num}.map ${testfile}_${split_num}.ped time_${threshold}.best > new_time_${threshold}.best
                                                echo "861"
                                                mv new_time_${threshold}.best time_${threshold}.best
                                                awk -v a="$threshold" 'NR==FNR{c[FNR]=$13;beta_add[FNR]=$6;type[FNR]=$1;beta_type[FNR]=$2;point[FNR]=$5;rs[FNR]=$10}ARGIND==2{for(snp=1;snp<=1000;snp++){if($2==rs[snp]){snp_site[FNR]=snp}}}ARGIND==3{human=FNR;for(snp=1;snp<=1000;snp++){if(point[snp_site[snp]]<a){if($(2*snp+5)==c[snp_site[snp]]&&$(2*snp+6)==c[snp_site[snp]]&&$(2*snp+5)!=0){sum[FNR]+=beta_add[snp_site[snp]]*2}else if($(2*snp+5)!=$(2*snp+6)||$(2*snp+5)==0){sum[FNR]+=beta_add[snp_site[snp]]}else{sum[FNR]+=0}}else if(type[snp_site[snp]]=="dom"){if(($(2*snp+5)==c[snp_site[snp]]&&$(2*snp+6)==c[snp_site[snp]]&&$(2*snp+5)!=0)||($(2*snp+5)!=$(2*snp+6))){sum[FNR]+=beta_type[snp_site[snp]]*2}else if($(2*snp+5)==0){sum[FNR]+=beta_type[snp_site[snp]]}else{sum[FNR]+=0}}else{if(($(2*snp+5)==c[snp_site[snp]]&&$(2*snp+6)==c[snp_site[snp]]&&$(2*snp+5)!=0)){sum[FNR]-=beta_type[snp_site[snp]]*2}else if($(2*snp+5)==0){sum[FNR]-=beta_type[snp_site[snp]]}else{sum[FNR]+=0}}}}ARGIND==4{print sum[FNR]+$1}' ${inputfile}.summary_and_bim_${split_num} ${validationfile}_${split_num}.map ${validationfile}_${split_num}.ped time_${threshold}_validation.best > new_time_${threshold}_validation.best
                                                mv new_time_${threshold}_validation.best time_${threshold}_validation.best
                                        } &
                                done
                                wait
                                if [ $with_origin = "True" ];then
                                        awk 'NR==FNR{c[FNR]=$13;beta_add[FNR]=$2;type_add[FNR]=$1;t_add[FNR]=$3;rs[FNR]=$10}ARGIND==2{for(snp=1;snp<=1000;snp++){if($2==rs[snp]){snp_site[FNR]=snp}}}ARGIND==3{human=FNR;for(snp=1;snp<=1000;snp++){if($(2*snp+5)==c[snp_site[snp]]&&$(2*snp+6)==c[snp_site[snp]]&&$(2*snp+5)!=0){sum[FNR]+=beta_add[snp_site[snp]]*2}else if($(2*snp+5)!=$(2*snp+6)||$(2*snp+5)==0){sum[FNR]+=beta_add[snp_site[snp]]}else{sum[FNR]+=0}}}ARGIND==4{print sum[FNR]+$1}' ${inputfile}.summary_and_bim_${split_num} ${testfile}_${split_num}.map ${testfile}_${split_num}.ped time.best > new_time.best
                                        mv new_time.best time.best
                                        awk 'NR==FNR{c[FNR]=$13;beta_add[FNR]=$2;type_add[FNR]=$1;t_add[FNR]=$3;rs[FNR]=$10}ARGIND==2{for(snp=1;snp<=1000;snp++){if($2==rs[snp]){snp_site[FNR]=snp}}}ARGIND==3{human=FNR;for(snp=1;snp<=1000;snp++){if($(2*snp+5)==c[snp_site[snp]]&&$(2*snp+6)==c[snp_site[snp]]&&$(2*snp+5)!=0){sum[FNR]+=beta_add[snp_site[snp]]*2}else if($(2*snp+5)!=$(2*snp+6)||$(2*snp+5)==0){sum[FNR]+=beta_add[snp_site[snp]]}else{sum[FNR]+=0}}}ARGIND==4{print sum[FNR]+$1}' ${inputfile}.summary_and_bim_${split_num} ${validationfile}_${split_num}.map ${validationfile}_${split_num}.ped time_validation.best > new_time_validation.best
                                        mv new_time_validation.best time_validation.best
                                fi
                        else
                                for threshold in $(seq $begin_site $step $end_site)
                                do
                                        {
                                                echo "872"
                                                awk -v a="$threshold" 'function max(a,b){if(a>b){return a}else{return b}};function min(a,b){if(a<b){return b}else{return a}};NR==FNR{c[FNR]=$13;beta_add[FNR]=$6;type[FNR]=$1;beta_type[FNR]=$2;point[FNR]=$5;rs[FNR]=$10}ARGIND==2{for(snp=1;snp<=1000;snp++){if($2==rs[snp]){snp_site[FNR]=snp}}}ARGIND==3{human=FNR;for(snp=1;snp<=1000;snp++){if($(2*snp+5)==c[snp_site[snp]]&&$(2*snp+6)==c[snp_site[snp]]&&$(2*snp+5)!=0){sum[FNR]+=beta_add[snp_site[snp]]*2*max(0,a-point[snp_site[snp]])}else if($(2*snp+5)!=$(2*snp+6)||$(2*snp+5)==0){sum[FNR]+=beta_add[snp_site[snp]]*max(0,a-point[snp_site[snp]])}else{sum[FNR]+=0}if(type[snp_site[snp]]=="dom"){if(($(2*snp+5)==c[snp_site[snp]]&&$(2*snp+6)==c[snp_site[snp]]&&$(2*snp+5)!=0)||($(2*snp+5)!=$(2*snp+6))){sum[FNR]+=beta_type[snp_site[snp]]*2*min(a,point[snp_site[snp]])}else if($(2*snp+5)==0){sum[FNR]+=beta_type[snp_site[snp]]*min(a,point[snp_site[snp]])}else{sum[FNR]+=0}}else{if(($(2*snp+5)==c[snp_site[snp]]&&$(2*snp+6)==c[snp_site[snp]]&&$(2*snp+5)!=0)){sum[FNR]-=beta_type[snp_site[snp]]*2*min(a,point[snp_site[snp]])}else if($(2*snp+5)==0){sum[FNR]-=beta_type[snp_site[snp]]*min(a,point[snp_site[snp]])}else{sum[FNR]+=0}}}}ARGIND==4{print sum[FNR]+$1}' ${inputfile}.summary_and_bim_${split_num} ${testfile}_${split_num}.map ${testfile}_${split_num}.ped time_${threshold}.best > new_time_${threshold}.best
                                                mv new_time_${threshold}.best time_${threshold}.best
                                                awk -v a="$threshold" 'function max(a,b){if(a>b){return a}else{return b}};function min(a,b){if(a<b){return b}else{return a}};NR==FNR{c[FNR]=$13;beta_add[FNR]=$6;type[FNR]=$1;beta_type[FNR]=$2;point[FNR]=$5;rs[FNR]=$10}ARGIND==2{for(snp=1;snp<=1000;snp++){if($2==rs[snp]){snp_site[FNR]=snp}}}ARGIND==3{human=FNR;for(snp=1;snp<=1000;snp++){if($(2*snp+5)==c[snp_site[snp]]&&$(2*snp+6)==c[snp_site[snp]]&&$(2*snp+5)!=0){sum[FNR]+=beta_add[snp_site[snp]]*2*max(0,a-point[snp_site[snp]])}else if($(2*snp+5)!=$(2*snp+6)||$(2*snp+5)==0){sum[FNR]+=beta_add[snp_site[snp]]*max(0,a-point[snp_site[snp]])}else{sum[FNR]+=0}if(type[snp_site[snp]]=="dom"){if(($(2*snp+5)==c[snp_site[snp]]&&$(2*snp+6)==c[snp_site[snp]]&&$(2*snp+5)!=0)||($(2*snp+5)!=$(2*snp+6))){sum[FNR]+=beta_type[snp_site[snp]]*2*min(a,point[snp_site[snp]])}else if($(2*snp+5)==0){sum[FNR]+=beta_type[snp_site[snp]]*min(a,point[snp_site[snp]])}else{sum[FNR]+=0}}else{if(($(2*snp+5)==c[snp_site[snp]]&&$(2*snp+6)==c[snp_site[snp]]&&$(2*snp+5)!=0)){sum[FNR]-=beta_type[snp_site[snp]]*2*min(a,point[snp_site[snp]])}else if($(2*snp+5)==0){sum[FNR]-=beta_type[snp_site[snp]]*min(a,point[snp_site[snp]])}else{sum[FNR]+=0}}}}ARGIND==4{print sum[FNR]+$1}' ${inputfile}.summary_and_bim_${split_num} ${validationfile}_${split_num}.map ${validationfile}_${split_num}.ped time_${threshold}_validation.best > new_time_${threshold}_validation.best
                                                mv new_time_${threshold}_validation.best time_${threshold}_validation.best
                                        } &
                                done
                                wait
                                if [ $with_origin = "True" ];then
                                        awk 'NR==FNR{c[FNR]=$13;beta_add[FNR]=$2;type_add[FNR]=$1;t_add[FNR]=$3;rs[FNR]=$10}ARGIND==2{for(snp=1;snp<=1000;snp++){if($2==rs[snp]){snp_site[FNR]=snp}}}ARGIND==3{human=FNR;for(snp=1;snp<=1000;snp++){if($(2*snp+5)==c[snp_site[snp]]&&$(2*snp+6)==c[snp_site[snp]]&&$(2*snp+5)!=0){sum[FNR]+=beta_add[snp_site[snp]]*2}else if($(2*snp+5)!=$(2*snp+6)||$(2*snp+5)==0){sum[FNR]+=beta_add[snp_site[snp]]}else{sum[FNR]+=0}}}ARGIND==4{print sum[FNR]+$1}' ${inputfile}.summary_and_bim_${split_num} ${testfile}_${split_num}.map ${testfile}_${split_num}.ped time.best > new_time.best
                                        mv new_time.best time.best
                                        awk 'NR==FNR{c[FNR]=$13;beta_add[FNR]=$2;type_add[FNR]=$1;t_add[FNR]=$3;rs[FNR]=$10}ARGIND==2{for(snp=1;snp<=1000;snp++){if($2==rs[snp]){snp_site[FNR]=snp}}}ARGIND==3{human=FNR;for(snp=1;snp<=1000;snp++){if($(2*snp+5)==c[snp_site[snp]]&&$(2*snp+6)==c[snp_site[snp]]&&$(2*snp+5)!=0){sum[FNR]+=beta_add[snp_site[snp]]*2}else if($(2*snp+5)!=$(2*snp+6)||$(2*snp+5)==0){sum[FNR]+=beta_add[snp_site[snp]]}else{sum[FNR]+=0}}}ARGIND==4{print sum[FNR]+$1}' ${inputfile}.summary_and_bim_${split_num} ${validationfile}_${split_num}.map ${validationfile}_${split_num}.ped time_validation.best > new_time_validation.best
                                        mv new_time_validation.best time_validation.best
                                fi
                        fi
                fi
                echo "882"
                symbol=0
                #rm ${inputfile}.summary_and_bim_${split_num} ${validationfile}_${split_num}.ped ${validationfile}_${split_num}.map ${validationfile}_${split_num}.bed ${validationfile}_${split_num}.bim ${validationfile}_${split_num}.fam ${num}_${split_num}_${threshold}_validation.best ${testfile}_${split_num}.ped ${testfile}_${split_num}.map ${testfile}_${split_num}.bed ${testfile}_${split_num}.bim ${testfile}_${split_num}.fam ${num}_${split_num}.best time_${threshold}.spearman time_${threshold}_validation.spearman
                rm ${validationfile}_${split_num}.ped ${validationfile}_${split_num}.map ${validationfile}_${split_num}.bed ${validationfile}_${split_num}.bim ${validationfile}_${split_num}.fam ${num}_${split_num}_${threshold}_validation.best ${testfile}_${split_num}.ped ${testfile}_${split_num}.map ${testfile}_${split_num}.bed ${testfile}_${split_num}.bim ${testfile}_${split_num}.fam ${num}_${split_num}.best time_${threshold}.spearman time_${threshold}_validation.spearman time_${threshold}.pearson time_${threshold}_validation.pearson time_${threshold}.r2 time_${threshold}_validation.r2 time_${threshold}.bce time_${threshold}_validation.bce
        done
else
        awk '$1~/add/{print $5}' ${testfile}.summary_and_bim > ${testfile}_add.txt
        awk '$1~/het/{print $5}' ${testfile}.summary_and_bim > ${testfile}_het.txt
        awk '$1~/rec/{print $5}' ${testfile}.summary_and_bim > ${testfile}_rec.txt
        awk '$1~/dom/{print $5}' ${testfile}.summary_and_bim > ${testfile}_dom.txt
        plink --bfile ${inputfile} --extract ${testfile}_add.txt --make-bed --out ${testfile}_add
        plink --bfile ${inputfile} --extract ${testfile}_het.txt --make-bed --out ${testfile}_het
        plink --bfile ${inputfile} --extract ${testfile}_rec.txt --make-bed --out ${testfile}_rec
        plink --bfile ${inputfile} --extract ${testfile}_dom.txt --make-bed --out ${testfile}_dom
        if [ $logistic = "True" ] ;then
                for threshold in $(seq 1 180)
                do
                        p_threshold=$(echo $input | awk 'NR=="'$threshold'"{print $0}' lower.txt)
                        /home/zhenghoufengLab/gaisirui/R-4.0.4/bin/Rscript	/home/zhenghoufengLab/gaisirui/PRSice.R	--dir	/home/zhenghoufengLab/gaisirui/R-4.0.4	--prsice	/home/zhenghoufengLab/gaisirui/PRSice_linux	--base	./${inputfile}.assoc.logistic	--target ${testfile}_dom	--pheno	${inputfile}.pheno	--thread	1	--beta					--binary-target	F	--stat	OR	--cov	${inputfile}.sex	--cov-col	@cov1	--cov-factor	@cov1	--all-score -i 1 -l $p_threshold --model dom --score std --out ${inputfile}_dom
                        /home/zhenghoufengLab/gaisirui/R-4.0.4/bin/Rscript	/home/zhenghoufengLab/gaisirui/PRSice.R	--dir	/home/zhenghoufengLab/gaisirui/R-4.0.4	--prsice	/home/zhenghoufengLab/gaisirui/PRSice_linux	--base	./${inputfile}.assoc.logistic	--target ${testfile}_add	--pheno	${inputfile}.pheno	--thread	1	--beta					--binary-target	F	--stat	OR	--cov	${inputfile}.sex	--cov-col	@cov1	--cov-factor	@cov1	--all-score -i 1 -l $p_threshold --model add --score std --out ${inputfile}_add
                        /home/zhenghoufengLab/gaisirui/R-4.0.4/bin/Rscript	/home/zhenghoufengLab/gaisirui/PRSice.R	--dir	/home/zhenghoufengLab/gaisirui/R-4.0.4	--prsice	/home/zhenghoufengLab/gaisirui/PRSice_linux	--base	./${inputfile}.assoc.logistic	--target ${testfile}_rec	--pheno	${inputfile}.pheno	--thread	1	--beta					--binary-target	F	--stat	OR	--cov	${inputfile}.sex	--cov-col	@cov1	--cov-factor	@cov1	--all-score -i 1 -l $p_threshold --model rec --score std --out ${inputfile}_rec
                        /home/zhenghoufengLab/gaisirui/R-4.0.4/bin/Rscript	/home/zhenghoufengLab/gaisirui/PRSice.R	--dir	/home/zhenghoufengLab/gaisirui/R-4.0.4	--prsice	/home/zhenghoufengLab/gaisirui/PRSice_linux	--base	./${inputfile}.assoc.logistic	--target ${testfile}_het	--pheno	${inputfile}.pheno	--thread	1	--beta					--binary-target	F	--stat	OR	--cov	${inputfile}.sex	--cov-col	@cov1	--cov-factor	@cov1	--all-score -i 1 -l $p_threshold --model het --score std --out ${inputfile}_het
                        paste ${inputfile}_dom.best ${inputfile}_add.best ${inputfile}_rec.best ${inputfile}_het.best > ${inputfile}_all_${threshold}.best
                        awk '{print $1"\t"$4+$8+$12+$16}' ${inputfile}_all_${threshold}.best ${inputfile}_sum.best
                        paste ${inputfile}_sum.best ${testfile}_${split_num}.fam > time.spearman
                        if [ $logistic = "False" ];then
                                relation=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import r2;print(r2.ReadData("time.spearman"))'`
                                echo "the p_value trying is"$threshold
                                echo "the p_value trying is"$threshold  >> $time_point.log
                                if [ $(echo "$old_relation <= $relation"|bc) = 1  ];then
                                        old_relation=$relation
                                        echo "the "$threshold" threshold got a better result as "$old_relation
                                        echo "the "$threshold" threshold got a better result as "$old_relation   >> $time_point.log
                                else
                                        echo "the "$threshold" threshold did not get a better result"
                                        echo "the "$threshold" threshold did not get a better result"   >> $time_point.log
                                fi
                        else
                                relation=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import bce;print(bce.ReadData("time.spearman"))'`
                                echo "the snp trying is"$threshold
                                echo "the snp trying is"$threshold    >> $time_point.log
                                if [ $(echo "$old_relation <= $relation"|bc) = 1  ];then
                                        old_relation=$relation
                                        echo "the "$threshold" threshold got a better result as "$old_relation
                                        echo "the "$threshold" threshold got a better result as "$old_relation     >> $time_point.log
                                        cp ${inputfile}_sum.best ${output_file}.prs
                                else
                                        echo "the "$threshold" snp did not get a better result"
                                        echo "the "$threshold" snp did not get a better result"     >> $time_point.log
                                fi
                        fi
                        rm ${num}_${split_num}.best ${num}_${split_num}.txt ${num}_${split_num}.ped ${num}_${split_num}.map ${num}_${split_num}.fam ${num}_${split_num}.bim ${num}_${split_num}.bed
                done
        else
                for threshold in $(seq 1 180)
                do
                        p_threshold=$(echo $input | awk 'NR=="'$threshold'"{print $0}' lower.txt)
                        /home/zhenghoufengLab/gaisirui/R-4.0.4/bin/Rscript	/home/zhenghoufengLab/gaisirui/PRSice.R	--dir	/home/zhenghoufengLab/gaisirui/R-4.0.4	--prsice	/home/zhenghoufengLab/gaisirui/PRSice_linux	--base	./${inputfile}.assoc.linear	--target ${testfile}_dom	--pheno	${inputfile}.pheno	--thread	1	--beta					--binary-target	F	--stat	BETA	--cov	${inputfile}.sex	--cov-col	@cov1	--cov-factor	@cov1	--all-score -i 1 -l $p_threshold --model dom --score std --out ${inputfile}_dom
                        /home/zhenghoufengLab/gaisirui/R-4.0.4/bin/Rscript	/home/zhenghoufengLab/gaisirui/PRSice.R	--dir	/home/zhenghoufengLab/gaisirui/R-4.0.4	--prsice	/home/zhenghoufengLab/gaisirui/PRSice_linux	--base	./${inputfile}.assoc.linear	--target ${testfile}_add	--pheno	${inputfile}.pheno	--thread	1	--beta					--binary-target	F	--stat	BETA	--cov	${inputfile}.sex	--cov-col	@cov1	--cov-factor	@cov1	--all-score -i 1 -l $p_threshold --model add --score std --out ${inputfile}_add
                        /home/zhenghoufengLab/gaisirui/R-4.0.4/bin/Rscript	/home/zhenghoufengLab/gaisirui/PRSice.R	--dir	/home/zhenghoufengLab/gaisirui/R-4.0.4	--prsice	/home/zhenghoufengLab/gaisirui/PRSice_linux	--base	./${inputfile}.assoc.linear	--target ${testfile}_rec	--pheno	${inputfile}.pheno	--thread	1	--beta					--binary-target	F	--stat	BETA	--cov	${inputfile}.sex	--cov-col	@cov1	--cov-factor	@cov1	--all-score -i 1 -l $p_threshold --model rec --score std --out ${inputfile}_rec
                        /home/zhenghoufengLab/gaisirui/R-4.0.4/bin/Rscript	/home/zhenghoufengLab/gaisirui/PRSice.R	--dir	/home/zhenghoufengLab/gaisirui/R-4.0.4	--prsice	/home/zhenghoufengLab/gaisirui/PRSice_linux	--base	./${inputfile}.assoc.linear	--target ${testfile}_het	--pheno	${inputfile}.pheno	--thread	1	--beta					--binary-target	F	--stat	BETA	--cov	${inputfile}.sex	--cov-col	@cov1	--cov-factor	@cov1	--all-score -i 1 -l $p_threshold --model het --score std --out ${inputfile}_het
                        paste ${inputfile}_dom.best ${inputfile}_add.best ${inputfile}_rec.best ${inputfile}_het.best > ${inputfile}_all.best
                        awk '{print $1"\t"$4+$8+$12+$16}' ${inputfile}_all.best ${inputfile}_sum.best
                        paste ${inputfile}_sum.best ${testfile}_${split_num}.fam > time.spearman
                        if [ $logistic = "False" ];then
                                relation=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import r2;print(r2.ReadData("time.spearman"))'`
                                echo "the p_value trying is"$threshold
                                echo "the p_value trying is"$threshold  >> $time_point.log
                                if [ $(echo "$old_relation <= $relation"|bc) = 1  ];then
                                        old_relation=$relation
                                        echo "the "$threshold" threshold got a better result as "$old_relation
                                        echo "the "$threshold" threshold got a better result as "$old_relation   >> $time_point.log
                                else
                                        echo "the "$threshold" threshold did not get a better result"
                                        echo "the "$threshold" threshold did not get a better result"   >> $time_point.log
                                fi
                        else
                                relation=`/home/zhenghoufengLab/gaisirui/miniconda3/bin/python -c 'import bce;print(bce.ReadData("time.spearman"))'`
                                echo "the snp trying is"$threshold
                                echo "the snp trying is"$threshold    >> $time_point.log
                                if [ $(echo "$old_relation <= $relation"|bc) = 1  ];then
                                        old_relation=$relation
                                        echo "the "$threshold" threshold got a better result as "$old_relation
                                        echo "the "$threshold" threshold got a better result as "$old_relation     >> $time_point.log
                                        cp ${inputfile}_sum.best ${output_file}.prs
                                else
                                        echo "the "$threshold" snp did not get a better result"
                                        echo "the "$threshold" snp did not get a better result"     >> $time_point.log
                                fi
                        fi
                        rm ${num}_${split_num}.best ${num}_${split_num}.txt ${num}_${split_num}.ped ${num}_${split_num}.map ${num}_${split_num}.fam ${num}_${split_num}.bim ${num}_${split_num}.bed
                done
        fi
fi
