#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2017/9/7 14:36
# @Author  : 
# @Site    : 
# @File    : cGcC_score.py
# @Software: PyCharmUPGIVEUP

#####################################
# cGcC[5]是一种新的打分系统，可以与标准模式结合从而更好地预测RNA中G-四联体的形成。
# cGcC分数分为两个部分，cG 分数是根据字符串s的长度n决定的： 
# 其中s表示序列字符串，Gs是s所有非空的连续的鸟嘌呤所组成的子串的集合，Gs(i)表示长度为i的子串的集合，|Gs(i)|表示集合Gs(i)的大小。
# 注：Gs中可以包含序列完全相同但位置不同的子串，例如：序列s=‘GGGCGGG’有两个’GGG’子串，所以|Gs(3)|=2。
# cG分数就定义为：cGs=∑_(i=1)^n▒(|Gs(i)|*10*i) 。
# cC分数定义相似为cCs=∑_(i=1)^n▒(|Cs(i)|*10*i) 。
# 也就是说，对于一段PG4序列，每一个G或是C给予分值10，每一段GG或是CC给予分值20，每一段GGG或是CCC给予分值30，以此类推。
# 例如序列s=’GGG’时，cGscore=3(G)*10+2(GG)*20+1(GGG)*30=100。s=’GG’时，cGscore=2(G)*10+1(GG)*20=40。
# 最后cGcC=(cG score)/(cC score)。
######################################


def lentoscore(len):
	'''
	_function:计算一段连续G子串或C子串的得分
	_输入：子串的长度__int型
	_输出：子串的分值__int型
	'''
    num = 1
    ite = len
    score = 0
    while ite > 0:
        score += ite * 10 * num
        ite = ite - 1
        num = num + 1
    return score


def cGcC(seq):
	'''
	_function:计算一段DNA序列或RNA序列的得分
	_输入：DNA序列或RNA序列__str型
	_输出：序列的cGcC得分__float型
	'''
    import re
    Gstandard = re.compile(r'G+')
    iterator = Gstandard.finditer(seq)
    cG = 0
    for i in iterator:
        cG += lentoscore(len(i.group()))
	## cG得分
    Cstandard = re.compile(r'C+')
    iterator = Cstandard.finditer(seq)
    cC = 0
    for i in iterator:
        cC += lentoscore(len(i.group()))
	## cC得分
    try:
        return cG / cC
    except ZeroDivisionError:
        return cG


if __name__ == '__main__':
    seq = 'ACAGGGASCGGGGCCCGGGAGVA'
    print(cGcC(seq))
