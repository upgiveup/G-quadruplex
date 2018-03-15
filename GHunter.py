#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2017/11/08 9:03
# @Author  : 
# @Site    : 
# @File    : GHunter.py
# @Software: PyCharmUPGIVEUP

##############################################
# G4Hunter是一种最新提出的G4打分方法，对序列进行打分时将序列分为两类，短序列或是长序列。
# 分类标准默认为25bp，也可人为设定。
# 在长序列时同时还需要设定G4Hx，一般设定在1-2之间，默认为1.5.
# 具体G4Hunter的分类原理可参考文档：G-四联体识别工具介绍和对比。
##############################################


def base_score(seq):
	'''
	_function：求G4Hunter分值的基础，对每一个碱基进行打分
	_input：DNA或RNA序列_str型
	_output：序列中每一个碱基的分值组成的列表_list型
	'''
    import re
    Gstandard = re.compile(r'G+')
    iterator_G = Gstandard.finditer(seq)
    Cstandard = re.compile(r'C+')
    iterator_C = Cstandard.finditer(seq)
    c = [0] * len(seq)
    for i in iterator_G:
        for j in range(i.span()[0], i.span()[1]):
            c[j] = min(len(i.group()), 4)
    for i in iterator_C:
        for j in range(i.span()[0], i.span()[1]):
            c[j] = -min(len(i.group()), 4)
    return c


def G_Hunter_short(seq):
    '''
	_function：求G4Hunter短序列的分值
	_input：DNA或RNA序列_str型
	_output：G4Hunter短序列的分值_float型
	'''
    score = base_score(seq)
    return sum(score) / len(seq)


def G_Hunter_long(seq, k, threshold):
	'''
	_function：求G4Hunter长序列的分值
	_input：seq: DNA或RNA序列_str型
			k：长短序列的分类标准_int型
			threshold：G4Hx_float_型
	_output：G4Hunter长序列的分值_float型
	'''
    score_base = base_score(seq)
    score = []
    for i in range(len(score_base) - k + 1):
        score.append(sum(score_base[i:i + k]) / k)
    loc = []
    for i in range(len(score)):
        if score[i] >= threshold:
            loc.append(i)
    try:
        tmp = loc[0] + k - 1
    except IndexError:
        print('No sequence meets the criterion!')
        return [], [], []
    sloc = []
    eloc = []
    sloc.append(loc[0])
    for i in loc:
        if i <= tmp:
            tmp = i + k - 1
        else:
            eloc.append(tmp)
            sloc.append(i)
            tmp = i + k - 1
    eloc.append(tmp)
    return score_base, sloc, eloc


def fusion(score, sloc, eloc):
	'''
	_function：将长序列中找到的序列精简为要输出的G4序列
	_input：score: 长序列中每一个碱基的分值_list型
			sloc：要精简的序列开始位点_int型
			eloc：要精简的序列结束位点_int型
	_output：sloc: 精简后的序列开始位点_int型
			 eloc: 精简后的序列结束位点_int型
	'''
    for i in range(len(sloc)):
        while score[sloc[i]] <= 0:
            sloc[i] += 1
    for i in range(len(eloc)):
        while score[eloc[i]] <= 0:
            eloc[i] -= 1
    return sloc, eloc


def G_Hunter(seq, k=25, threshold=1.5):
	'''
	_function：求G4Hunter序列的分值
	_input：seq: DNA或RNA序列_str型
			k：长短序列的分类标准_int型
			threshold：G4Hx_float_型
	_output：G4Hunter短序列的分值_float型
			G4Hunter长序列中的G4序列子串及其分值_dict型
	'''
    if len(seq) <= k or k == -1:
        return G_Hunter_short(seq)
    else:
        score, sloc, eloc = G_Hunter_long(seq, k, threshold)
        sloc, eloc = fusion(score, sloc, eloc)
        G4 = []
        score_dict = {}
        index_G4=0
        for i in range(len(sloc)):
            G4.append(seq[sloc[i]:(eloc[i] + 1)])
            scor=sum(score[sloc[i]:(eloc[i] + 1)])/(eloc[i] + 1-sloc[i])
            score_dict[G4[index_G4]]=scor
            index_G4+=1
        return score_dict


if __name__ == '__main__':
    seq = 'GAGGACAAGGAGGTGCGAGGAAAGGGGTTGGGGGATGGTCCCACAGGCAGCCACACCTGAGGCGTGGGCGGCCGGTAGGAGCTGGGGGAGGGCG' \
          'GGGAGAAGAGGGGTTTCTGTGTGTAG'
    #seq = 'GGGCCCCCCGGGCCCCCCGGGCCCCCCGGG'
    print(len(seq))
    print(G_Hunter(seq, k=35, threshold=0))
