本仓库包含一个用于构造一对 **正交拉丁方 (Orthogonal Latin Squares)** 的 C++ 实现。

一对正交拉丁方是满足以下条件的两个拉丁方：

> 若 $A$ 和 $B$ 为阶数 $n$ 的拉丁方，则对于任意 $(i,j)\ne(k,l)$，有 $(A_{ij},B_{ij})\ne(A_{kl},B_{kl})$。

本项目基于[这里](https://www.zhihu.com/question/1964254729556170187/answer/1964473533217347184)中给出的构造方法，实现了其构造。

代码编写者：chenyuxuan2009

This repository contains a C++ implementation for constructing a pair of **Orthogonal Latin Squares (OLS)**.

A pair of orthogonal Latin squares satisfies the following condition:

> If $A$ and $B$ are Latin squares of order $n$, then for any $(i, j) \ne (k, l)$, we have $(A_{ij}, B_{ij}) \ne (A_{kl}, B_{kl})$.

This project implements the construction method described [here](https://www.zhihu.com/question/1964254729556170187/answer/1964473533217347184).

Code author: chenyuxuan2009

