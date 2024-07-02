nodes=("in006" "in007" "in008")
#
#"in005" "in006" "in007" "in008" "in009" "in010"
node_count=${#nodes[@]}
counter=0
# 循环生成.sh文件
for j in {4..14}
do
echo $j
    # 获取当前节点
    node="${nodes[$counter % $node_count]}"
    # 定义.sh文件名
    filename="script_${j}.sh"
    # 创建.sh文件并写入内容
    cat <<EOF > "$filename"
#!/bin/bash
#SBATCH -e scPagwas_$j.err
#SBATCH -o scPagwas_$j.out
#SBATCH -J scPagwas_$j
#SBATCH -w $node
#SBATCH --mem=200000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate myr
Rscript /share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/1.r $j
EOF
    # 增加节点计数器
    counter=$((counter + 1))
    sbatch $filename
done
#$node

nodes=("in006" "in007" "in008")
#
#"in005" "in006" "in007" "in008" "in009" "in010"
node_count=${#nodes[@]}
counter=0
# 循环生成.sh文件
for j in {1..11}
do
echo $j
    # 获取当前节点
    node="${nodes[$counter % $node_count]}"
    # 定义.sh文件名
    filename="script_${j}.sh"
    # 创建.sh文件并写入内容
    cat <<EOF > "$filename"
#!/bin/bash
#SBATCH -e scPagwas_$j.err
#SBATCH -o scPagwas_$j.out
#SBATCH -J scPagwas_$j
#SBATCH -w $node
#SBATCH --mem=200000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate myr
Rscript /share/pub/dengcy/Cancer_Gwas/src1.0/4.0.Run_scPagwas/1.r $j
EOF
    # 增加节点计数器
    counter=$((counter + 1))
    sbatch $filename
done