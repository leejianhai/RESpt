import subprocess
import sys
import os
from pathlib import Path

def check_dependencies():
    """检查并安装必要的Python包"""
    required_packages = ['pandas', 'biopython']
    
    for package in required_packages:
        try:
            __import__(package)
        except ImportError:
            print(f"正在安装 {package}...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])
    
    print("所有依赖包已安装完成")

# 首先检查并安装依赖
check_dependencies()

# 然后导入其他包
import pandas as pd
from Bio import SeqIO
import json
import argparse

class GenomeAnalyzer:
    def __init__(self, input_genome, outdir):
        self.input_genome = input_genome
        self.outdir = Path(outdir)
        
        # 设置默认数据库和输出文件路径
        self.setup_paths()
        
    def setup_paths(self):
        """设置各种文件路径"""
        # 创建输出目录
        self.outdir.mkdir(parents=True, exist_ok=True)
        
        # 数据库路径 - 使用绝对路径
        self.db_dir = Path('../RESpt')
        self.subject_fasta_ISTnIN = self.db_dir / 'IS_Tn_IN/tncentral_integrall_isfinder.fa'
        self.subject_fasta_RESfinder = self.db_dir / 'resfinder/resfinder.fna'
        
        # Platon数据库路径
        self.platon_db = Path('../RESpt/platon')
        
        # 输出文件路径
        self.output_file_ISTnIN = self.outdir / 'ISTnIN.tsv'
        self.output_file_RESfinder = self.outdir / 'RESfinder.tsv'
        self.prokka_dir = self.outdir / 'prokka_results'
        self.platon_dir = self.outdir / 'platon_results'
        self.phispy_dir = self.outdir / 'phispy_results'

    def check_conda_env(self, env_name):
        """检查conda环境是否存在"""
        result = subprocess.run(['conda', 'env', 'list'], capture_output=True, text=True)
        return env_name in result.stdout

    def check_sequence_length(self, fasta_file):
        """检查FASTA文件中的序列长度"""
        min_length = float('inf')
        for record in SeqIO.parse(fasta_file, "fasta"):
            min_length = min(min_length, len(record.seq))
        return min_length

    def check_output_exists(self, output_path, description):
        """检查输出文件是否存在
        
        参数:
            output_path: 输出文件路径
            description: 分析步骤描述
        
        返回:
            bool: 如果文件存在返回True，否则返回False
        """
        if output_path.exists():
            print(f"{description}结果已存在：{output_path}")
            return True
        return False

    def run_prokka(self):
        """运行Prokka注释"""
        if not self.check_conda_env('prokka'):
            print("错误：未找到prokka conda环境")
            return
            
        prefix = Path(self.input_genome).stem
        self.prokka_dir.mkdir(parents=True, exist_ok=True)
        
        # 检查GBK文件是否已存在
        prokka_gbk = self.prokka_dir / f"{prefix}.gbk"
        if self.check_output_exists(prokka_gbk, "Prokka"):
            return
        
        try:
            subprocess.run([
                'conda', 'run', '-n', 'prokka', 'prokka',
                str(self.input_genome),
                '--outdir', str(self.prokka_dir),
                '--prefix', prefix,
                '--force'
            ], check=True)
            print(f"Prokka结果已保存至：{self.prokka_dir}")
        except subprocess.CalledProcessError as e:
            print(f"运行Prokka时出错：{e}")

    def filter_blast_results(self, blast_result_file, query_fasta, subject_fasta, output_file):
        """过滤BLAST结果
        
        参数:
            blast_result_file: BLAST输出文件
            query_fasta: 查询序列文件
            subject_fasta: 目标序列文件
            output_file: 过滤后的输出文件
        """
        # 定义BLAST输出的列名
        columns = [
            'query_id', 'subject_id', 'identity', 'alignment_length',
            'mismatches', 'gap_opens', 'query_start', 'query_end',
            'subject_start', 'subject_end', 'evalue', 'bit_score'
        ]
        
        # 读取BLAST结果
        blast_df = pd.read_csv(blast_result_file, sep='\t', names=columns)
        
        # 获取查询序列长度
        query_lengths = {}
        for record in SeqIO.parse(query_fasta, "fasta"):
            query_lengths[record.id] = len(record.seq)
            
        # 获取目标序列长度
        subject_lengths = {}
        for record in SeqIO.parse(subject_fasta, "fasta"):
            subject_lengths[record.id] = len(record.seq)
            
        # 添加序列长度到DataFrame
        blast_df['query_length'] = blast_df['query_id'].map(query_lengths)
        blast_df['subject_length'] = blast_df['subject_id'].map(subject_lengths)
        
        # 计算覆盖率
        blast_df['query_coverage'] = (blast_df['alignment_length'] / blast_df['query_length']) * 100
        blast_df['subject_coverage'] = (blast_df['alignment_length'] / blast_df['subject_length']) * 100
        
        # 过滤结果 - 覆盖率和相似度
        filtered_df = blast_df[
            (blast_df['query_coverage'] >= 80) &
            (blast_df['subject_coverage'] >= 80) &  # 添加subject_coverage过滤
            (blast_df['identity'] >= 80)
        ]
        
        # 对每个query_id只保留identity最高的结果
        filtered_df = filtered_df.sort_values('identity', ascending=False)
        filtered_df = filtered_df.loc[filtered_df.groupby('query_id')['identity'].idxmax()]
        
        # 选择输出列
        output_columns = [
            'query_id', 'subject_id', 'identity', 'alignment_length',
            'query_start', 'query_end', 'query_coverage', 
            'subject_coverage',  # 添加subject_coverage到输出
            'evalue', 'bit_score'
        ]
        filtered_df = filtered_df[output_columns]
        
        # 按identity降序排序
        filtered_df = filtered_df.sort_values('identity', ascending=False)
        
        # 保存过滤后的结果
        filtered_df.to_csv(output_file, sep='\t', index=False)
        print(f"过滤后的BLAST结果已保存至：{output_file}")
        
        return filtered_df

    def run_blast(self, query_fasta, data_base, blast_result_file):
        """运行BLAST分析"""
        if not self.check_conda_env('blast'):
            print("错误：未找到blast conda环境")
            return
            
        # 检查过滤后的BLAST结果是否已存在
        filtered_output = Path(str(blast_result_file).replace('.tsv', '_filtered.tsv'))
        if self.check_output_exists(filtered_output, "BLAST"):
            return
            
        try:
            # 运行BLAST
            subprocess.run([
                'conda', 'run', '-n', 'blast', 'blastn',
                '-query', query_fasta,
                '-db', data_base,
                '-evalue', '1e-5',
                '-out', blast_result_file,
                '-outfmt', '6'
            ], check=True)
            print(f"BLAST结果已保存至：{blast_result_file}")
            
            # 过滤结果
            self.filter_blast_results(blast_result_file, query_fasta, data_base, filtered_output)
            
        except subprocess.CalledProcessError as e:
            print(f"运行BLAST时出错：{e}")

    def run_platon(self):
        """运行Platon分析"""
        if not self.check_conda_env('platon'):
            print("错误：未找到platon conda环境")
            return
            
        if not self.platon_db.exists():
            print(f"错误：未找到Platon数据库：{self.platon_db}")
            return
            
        prefix = Path(self.input_genome).stem
        self.platon_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            subprocess.run([
                'conda', 'run', '-n', 'platon', 'platon',
                str(self.input_genome),
                '--db', str(self.platon_db),
                '--output', str(self.platon_dir),
                '--prefix', prefix
            ], check=True)
            print(f"Platon结果已保存至：{self.platon_dir}")
            
            # 检查并解析JSON文件
            json_file = self.platon_dir / f"{prefix}.json"
            if not json_file.exists():
                print(f"警告：未找到Platon的JSON输出文件：{json_file}")
                # 检查其他可能的文件位置
                alt_json_file = Path(f"{prefix}.json")
                if alt_json_file.exists():
                    json_file = alt_json_file
                    print(f"在当前目录找到JSON文件：{json_file}")
            
            if json_file.exists():
                self.parse_platon_json(
                    json_file,
                    self.outdir / 'platon_results_summary.tsv'
                )
            else:
                print("无法找到Platon的JSON输出文件")
            
        except subprocess.CalledProcessError as e:
            print(f"运行Platon时出错：{e}")

    def run_phispy(self):
        """运行PhiSpy分析"""
        if not self.check_conda_env('phispy'):
            print("错误：未找到phispy conda环境")
            return
            
        prefix = Path(self.input_genome).stem
        self.phispy_dir.mkdir(parents=True, exist_ok=True)
        
        # 使用Prokka的gbk文件作为输入
        prokka_gbk = self.prokka_dir / f"{prefix}.gbk"
        if not prokka_gbk.exists():
            print(f"错误：未找到Prokka的GBK文件：{prokka_gbk}")
            return
            
        try:
            subprocess.run([
                'conda', 'run', '-n', 'phispy', 'phispy',
                str(prokka_gbk),
                '-o', str(self.phispy_dir),
                '-p', prefix,
                '--threads', '4'
            ], check=True)
            print(f"PhiSpy结果已保存至：{self.phispy_dir}")
            
            # 检查并处理输出文件
            prophage_files = list(self.phispy_dir.glob('*prophage_coordinates.tsv'))
            if prophage_files:
                prophage_file = prophage_files[0]  # 使用找到的第一个文件
                print(f"找到PhiSpy输出文件：{prophage_file}")
                
                # 读取并处理噬菌体预测结果
                try:
                    phage_df = pd.read_csv(prophage_file, sep='\t')
                    # 保存一个格式化的版本
                    output_file = self.outdir / 'phispy_results_summary.tsv'
                    phage_df.to_csv(output_file, sep='\t', index=False)
                    print(f"噬菌体预测结果已保存至：{output_file}")
                except Exception as e:
                    print(f"处理PhiSpy输出文件时出错：{e}")
            else:
                print("未找到任何PhiSpy输出文件")
            
        except subprocess.CalledProcessError as e:
            print(f"运行PhiSpy时出错：{e}")

    def parse_gff_file(self, gff_file):
        """解析GFF文件获取基因位置信息
        
        返回字典：
        {
            gene_id: {
                'contig': contig_id,
                'start': start_pos,
                'end': end_pos,
                'strand': strand
            }
        }
        """
        gene_locations = {}
        try:
            with open(gff_file) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    parts = line.strip().split('\t')
                    if len(parts) < 9:
                        continue
                        
                    if parts[2] == 'CDS':  # 只处理CDS条目
                        contig = parts[0]
                        start = int(parts[3])
                        end = int(parts[4])
                        strand = parts[6]
                        
                        # 从属性字段中提取基因ID
                        attributes = dict(item.split('=') for item in parts[8].split(';') if '=' in item)
                        if 'ID' in attributes:
                            gene_id = attributes['ID']
                            gene_locations[gene_id] = {
                                'contig': contig,
                                'start': start,
                                'end': end,
                                'strand': strand
                            }
            
            print(f"从GFF文件中提取了 {len(gene_locations)} 个基因的位置信息")
        except Exception as e:
            print(f"解析GFF文件时出错：{e}")
        
        return gene_locations

    def add_genome_locations(self, blast_file, gff_file, output_file, is_istnin=False):
        """添加基因组位置信息到BLAST结果"""
        # 读取BLAST结果
        blast_df = pd.read_csv(blast_file, sep='\t')
        
        # 从GFF文件中提取基因位置信息
        gene_locations = self.parse_gff_file(gff_file)
        
        # 添加位置信息到DataFrame
        def get_location(query_id):
            # 处理Prokka的基因ID格式（如果需要）
            gene_id = query_id
            if not gene_id in gene_locations:
                # 尝试添加'CDS_'前缀（Prokka的GFF文件中可能使用这种格式）
                gene_id = f"CDS_{query_id}"
            
            loc = gene_locations.get(gene_id, {})
            if not loc:
                print(f"Warning: No location found for {query_id}")
            return loc
        
        # 添加位置信息
        blast_df['contig'] = blast_df['query_id'].map(lambda x: get_location(x).get('contig'))
        blast_df['genome_start'] = blast_df['query_id'].map(lambda x: get_location(x).get('start'))
        blast_df['genome_end'] = blast_df['query_id'].map(lambda x: get_location(x).get('end'))
        blast_df['strand'] = blast_df['query_id'].map(lambda x: get_location(x).get('strand'))
        
        # 如果是ISTnIN结果，添加类型标记
        if is_istnin:
            def get_element_type(subject_id):
                if 'ISFinder' in subject_id:
                    return 'IS'
                elif 'INTEGRALL' in subject_id:
                    return 'IN'
                else:
                    return 'Tn'
            
            blast_df['element_type'] = blast_df['subject_id'].apply(get_element_type)
            
            # 更新列顺序以包含element_type
            columns = [
                'query_id', 'subject_id', 'identity', 'alignment_length',
                'query_start', 'query_end', 'query_coverage', 'subject_coverage',
                'contig', 'genome_start', 'genome_end', 'strand',
                'evalue', 'bit_score', 'element_type'
            ]
        else:
            columns = [
                'query_id', 'subject_id', 'identity', 'alignment_length',
                'query_start', 'query_end', 'query_coverage', 'subject_coverage',
                'contig', 'genome_start', 'genome_end', 'strand',
                'evalue', 'bit_score'
            ]
        
        result_df = blast_df[columns]
        
        # 检查是否成功添加了位置信息
        location_added = result_df['genome_start'].notna().any()
        if not location_added:
            print("警告：未能成功添加位置信息！")
            print("示例记录：")
            print(result_df.head())
            print("\n原始BLAST文件中的query_id示例：")
            print(blast_df['query_id'].head())
            print("\n提取到的基因位置信息例：")
            print(list(gene_locations.keys())[:5])
        
        # 保存结果
        result_df.to_csv(output_file, sep='\t', index=False)
        print(f"已添加基因组位置信息并保存至：{output_file}")
        
        return result_df

    def analyze_resistance_genes_context(self):
        """分析耐药基因的基因组环境"""
        try:
            # 读取所有必要的文件
            res_file = self.outdir / 'RESfinder_filtered_with_locations.tsv'
            istnin_file = self.outdir / 'ISTnIN_filtered_with_locations.tsv'
            phage_file = self.outdir / 'phispy_results_summary.tsv'
            plasmid_file = self.outdir / 'platon_results_summary.tsv'
            
            if not all(f.exists() for f in [res_file, istnin_file]):
                print("缺少必要的输入文件")
                return
            
            # 读取文件
            res_df = pd.read_csv(res_file, sep='\t')
            istnin_df = pd.read_csv(istnin_file, sep='\t')
            phage_df = pd.read_csv(phage_file, sep='\t') if phage_file.exists() else pd.DataFrame()
            plasmid_df = pd.read_csv(plasmid_file, sep='\t') if plasmid_file.exists() else pd.DataFrame()
            
            # 分离不同类型的移动元件
            is_df = istnin_df[istnin_df['element_type'] == 'IS']
            in_df = istnin_df[istnin_df['element_type'] == 'IN']
            tn_df = istnin_df[istnin_df['element_type'] == 'Tn']
            
            # 创建结果列表
            results = []
            
            # 分析每个耐药基因
            for _, res_gene in res_df.iterrows():
                result = {
                    'resistance_gene': res_gene['query_id'],
                    'resistance_gene_name': res_gene['subject_id'],
                    'contig': res_gene['contig'],
                    'start': res_gene['genome_start'],
                    'end': res_gene['genome_end'],
                    'identity': res_gene['identity'],
                    'on_plasmid': 'No',
                    'on_phage': 'No',
                    'in_integron': 'No',
                    'in_transposon': 'No',
                    'flanking_IS': 'No',
                    'flanking_IN': 'No',
                    'mobile_potential': 'No',
                    'mobile_evidence': []
                }
                
                # 检查是否在质粒上
                if not plasmid_df.empty:
                    for _, plasmid in plasmid_df.iterrows():
                        if (res_gene['contig'] == plasmid['contig_id'] and
                            res_gene['genome_start'] >= plasmid['hit_contig_start'] and
                            res_gene['genome_end'] <= plasmid['hit_contig_end']):
                            result['on_plasmid'] = 'Yes'
                            result['mobile_evidence'].append(f"Located on plasmid {plasmid['hit_plasmid_id']}")
                            break
                
                # 检查是否在噬菌体区域
                if not phage_df.empty:
                    for _, phage in phage_df.iterrows():
                        if (res_gene['contig'] == phage['Contig'] and
                            res_gene['genome_start'] >= phage['Start'] and
                            res_gene['genome_end'] <= phage['End']):
                            result['on_phage'] = 'Yes'
                            result['mobile_evidence'].append('Located in prophage region')
                            break
                # 检查是否在转座子范围内
                for _, tn in tn_df.iterrows():
                    if (res_gene['contig'] == tn['contig'] and
                        res_gene['genome_start'] >= tn['genome_start'] and
                        res_gene['genome_end'] <= tn['genome_end']):
                        result['in_transposon'] = 'Yes'
                        result['mobile_evidence'].append(f"Located within transposon {tn['subject_id']}")
                        break
                
                # 检查是否在整合子范围内
                for _, integron in in_df.iterrows():
                    if (res_gene['contig'] == integron['contig'] and
                        res_gene['genome_start'] >= integron['genome_start'] and
                        res_gene['genome_end'] <= integron['genome_end']):
                        result['in_integron'] = 'Yes'
                        result['mobile_evidence'].append(f"Located within integron {integron['subject_id']}")
                        break

                
                # 检查10kb范围内的移动元件
                window = 10000
                
                # 检查左右两侧是否有相同的IN
                # 检查上游10kb范围内的IN
                left_in = in_df[
                    (in_df['contig'] == res_gene['contig']) &
                    (in_df['genome_end'] <= res_gene['genome_start']) &
                    (in_df['genome_end'] >= res_gene['genome_start'] - window)
                ]
                
                # 检查下游10kb范围内的IN
                right_in = in_df[
                    (in_df['contig'] == res_gene['contig']) &
                    (in_df['genome_start'] >= res_gene['genome_end']) &
                    (in_df['genome_start'] <= res_gene['genome_end'] + window)
                ]
                
                # 获取上下游IN类型
                left_in_types = set(left_in['subject_id'])
                right_in_types = set(right_in['subject_id'])
                
                # 检查是否存在相同的IN
                common_in = left_in_types.intersection(right_in_types)
                if common_in:
                    result['flanking_IN'] = 'Yes'
                    result['mobile_evidence'].append(f"两侧存在相同的整合子: {';'.join(common_in)}")
                    result['mobile_potential'] = 'Yes'
                
                # 检查上游10kb范围内的IS
                left_is = is_df[
                    (is_df['contig'] == res_gene['contig']) &
                    (is_df['genome_end'] <= res_gene['genome_start']) &
                    (is_df['genome_end'] >= res_gene['genome_start'] - window)
                ]
                
                # 检查下游10kb范围内的IS
                right_is = is_df[
                    (is_df['contig'] == res_gene['contig']) &
                    (is_df['genome_start'] >= res_gene['genome_end']) &
                    (is_df['genome_start'] <= res_gene['genome_end'] + window)
                ]
                
                # 获取上下游IS类型
                left_is_types = set(left_is['subject_id'])
                right_is_types = set(right_is['subject_id'])
                
                # 检查是否存在相同的IS
                common_is = left_is_types.intersection(right_is_types)
                if common_is:
                    result['flanking_IS'] = 'Yes'
                    result['mobile_evidence'].append(f"两侧存在相同的IS序列: {';'.join(common_is)}")
                    result['mobile_potential'] = 'Yes'
                
                # 判断耐药基因转移可能性
                if (result['on_plasmid'] == 'Yes' or 
                    result['on_phage'] == 'Yes' or 
                    result['in_integron'] == 'Yes' or 
                    result['in_transposon'] == 'Yes' or 
                    result['flanking_IS'] == 'Yes' or 
                    result['flanking_IN'] == 'Yes'):
                    result['mobile_potential'] = 'Yes'
                
                results.append(result)
            
            # 转换为DataFrame
            results_df = pd.DataFrame(results)
            
            # 将列表转换为字符串
            results_df['mobile_evidence'] = results_df['mobile_evidence'].apply(lambda x: ' | '.join(x) if x else '')
            
            # 保存结果
            output_file = self.outdir / 'resistance_genes_context.tsv'
            results_df.to_csv(output_file, sep='\t', index=False)
            print(f"耐药基因环境分析结果已保存至：{output_file}")
            
            return results_df
            
        except Exception as e:
            print(f"分析耐药基因环境时出错：{e}")
            return None

    def run_analysis(self):
        """运行完整的分析流程"""
        print(f"开始分析基因组：{self.input_genome}")
        print(f"输出目录：{self.outdir}")
        
        # 首先运行Prokka，因为PhiSpy需要其输出
        self.run_prokka()
        
        # 运行BLAST分析
        prokka_ffn = self.prokka_dir / f"{Path(self.input_genome).stem}.ffn"
        prokka_gff = self.prokka_dir / f"{Path(self.input_genome).stem}.gff"
        
        if prokka_ffn.exists() and prokka_gff.exists():
            # 运行BLAST和过滤
            self.run_blast(str(prokka_ffn), str(self.subject_fasta_ISTnIN), str(self.output_file_ISTnIN))
            self.run_blast(str(prokka_ffn), str(self.subject_fasta_RESfinder), str(self.output_file_RESfinder))
            
            # 添加基因组位置信息
            istnin_filtered = Path(str(self.output_file_ISTnIN).replace('.tsv', '_filtered.tsv'))
            resfinder_filtered = Path(str(self.output_file_RESfinder).replace('.tsv', '_filtered.tsv'))
            
            if istnin_filtered.exists():
                self.add_genome_locations(
                    istnin_filtered,
                    prokka_gff,
                    self.outdir / 'ISTnIN_filtered_with_locations.tsv',
                    is_istnin=True  # 标记这是ISTnIN结果
                )
            
            if resfinder_filtered.exists():
                self.add_genome_locations(
                    resfinder_filtered,
                    prokka_gff,
                    self.outdir / 'RESfinder_filtered_with_locations.tsv',
                    is_istnin=False
                )
        else:
            print(f"错误：未找到Prokka输出文件")
        
        # 运行其他分析
        print("\n运行Platon分析...")
        self.run_platon()
        
        print("\n运行PhiSpy分析...")
        self.run_phispy()
        
        # 在所有分析完成后，添加耐药基因环境分析
        print("\n分析耐药基因环境...")
        self.analyze_resistance_genes_context()

    def parse_platon_json(self, json_file, output_file):
        """解析Platon的JSON文件并提取质粒信息
        
        参数:
            json_file: Platon生成的JSON文件路径
            output_file: 输出TSV文件路径
        """
        try:
            # 读取JSON文件
            with open(json_file) as f:
                data = json.load(f)
            
            # 存储所有提取的信息
            plasmid_info = []
            
            # 遍历每个contig
            for contig_id, contig_data in data.items():
                # 基本信息
                basic_info = {
                    'contig_id': contig_id,
                    'length': contig_data.get('length', ''),
                    'gc': contig_data.get('gc', ''),
                    'coverage': contig_data.get('coverage', ''),
                    'circular': contig_data.get('circular', False),
                    'plasmid_probability': contig_data.get('plasmid_probability', '')
                }
                
                # 提取质粒hits信息
                plasmid_hits = contig_data.get('plasmid_hits', [])
                for hit in plasmid_hits:
                    hit_info = basic_info.copy()
                    
                    # 添加hit特定信息
                    if 'plasmid' in hit:
                        hit_info.update({
                            'hit_plasmid_id': hit['plasmid'].get('id', ''),
                            'hit_plasmid_description': hit['plasmid'].get('description', ''),
                            'hit_plasmid_length': hit['plasmid'].get('length', ''),
                            'hit_identity': hit.get('identity', ''),
                            'hit_coverage': hit.get('coverage', ''),
                            'hit_contig_start': hit.get('contig_start', ''),
                            'hit_contig_end': hit.get('contig_end', ''),
                            'hit_plasmid_start': hit.get('plasmid_start', ''),
                            'hit_plasmid_end': hit.get('plasmid_end', '')
                        })
                    
                    plasmid_info.append(hit_info)
            
            # 如果没有找到任何质粒信息
            if not plasmid_info:
                print("未在JSON文件中找到质粒信息")
                return None
            
            # 转换为DataFrame
            df = pd.DataFrame(plasmid_info)
            
            # 设置列的顺序
            columns = [
                'contig_id', 'length', 'gc', 'coverage', 'circular', 
                'plasmid_probability', 'hit_plasmid_id', 'hit_plasmid_description',
                'hit_plasmid_length', 'hit_identity', 'hit_coverage',
                'hit_contig_start', 'hit_contig_end',
                'hit_plasmid_start', 'hit_plasmid_end'
            ]
            df = df[columns]
            
            # 保存为TSV文件
            df.to_csv(output_file, sep='\t', index=False)
            print(f"质粒信息已保存至：{output_file}")
            
            return df
            
        except Exception as e:
            print(f"解析Platon JSON文件时出错：{e}")
            return None

def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='基因组分析流程')
    parser.add_argument('input_genome', help='输入基因组文件路径')
    parser.add_argument('--outdir', default='../analyze', help='输出目录路径 (默认: ../analyze)')
    return parser.parse_args()

def main():
    args = parse_args()
    analyzer = GenomeAnalyzer(args.input_genome, args.outdir)
    analyzer.run_analysis()

if __name__ == "__main__":
    main()
 