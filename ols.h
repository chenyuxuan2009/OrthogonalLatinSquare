#ifndef PROJECT_OLS_H_
#define PROJECT_OLS_H_
#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <set>
#include <map>
#include <cstdint>
#include <cstring>
#include <climits>
#include <iostream>
namespace GF
{
    /*
    thanks to ChatGPT for implementing the GF(2^k) arithmetic
    */
    // 多项式乘法 mod f(x)，在 GF(2) 上
    uint64_t poly_mul_mod(uint64_t a, uint64_t b, uint64_t f, int k)
    {
        uint64_t res = 0;
        for (int i = 0; i < k; ++i)
            if ((b >> i) & 1)
                res ^= (a << i);
        for (int bit = (k << 1) - 2; bit >= k; --bit)
            if ((res >> bit) & 1)
                res ^= f << (bit - k);
        return res & ((1ULL << k) - 1);
    }

    // 多项式幂 mod f(x)：返回 (x^(2^p)) mod f
    uint64_t poly_square_mod(uint64_t x, uint64_t f, int k, int p)
    {
        uint64_t res = 2; // 初始为 x
        for (int i = 0; i < p; ++i)
            res = poly_mul_mod(res, res, f, k);
        return res;
    }

    // gcd over GF(2)[x]
    uint64_t poly_gcd(uint64_t a, uint64_t b)
    {
        while (b)
        {
            int da = 63 - __builtin_clzll(a);
            int db = 63 - __builtin_clzll(b);
            while (da >= db)
            {
                a ^= b << (da - db);
                da = 63 - __builtin_clzll(a);
            }
            std::swap(a, b);
        }
        return a;
    }

    // 判断多项式 f 是否不可约
    bool is_irreducible(uint64_t f, int k)
    {
        uint64_t x = 2; // x
        for (int i = 1; i <= k / 2; ++i)
        {
            x = poly_square_mod(x, f, k, 1);
            uint64_t g = poly_gcd(f, x ^ 2);
            if (g != 1)
                return false;
        }
        return true;
    }

    // 查表 + 自动求解
    uint64_t irreducible_poly(int k)
    {
        const std::map<int, uint64_t> table = {
            {1, 0b11},
            {2, 0b111},
            {3, 0b1011},
            {4, 0b10011},
            {5, 0b100101},
            {6, 0b1000011},
            {7, 0b10000011},
            {8, 0b100011101},
            {9, 0b1000010001},
            {10, 0b10000001001},
            {11, 0b100000000101},
            {12, 0b1000001010011},
            {13, 0b10000000011011},
            {14, 0b100010000000011},
            {15, 0b1100000000000001},
            {16, 0b10001000000001011}};
        auto it = table.find(k);
        if (it != table.end())
            return it->second;

        for (uint64_t f = (1ULL << k) | 1; f < (1ULL << (k + 1)); f += 2)
            if (is_irreducible(f, k))
                return f;
        throw std::runtime_error("irreducible polynomial not found");
    }

    uint64_t gf_mul(uint64_t a, uint64_t b, uint64_t poly, int k)
    {
        uint64_t res = 0;
        for (int i = 0; i < k; ++i)
            if ((b >> i) & 1)
                res ^= (a << i);
        for (int bit = (k << 1) - 2; bit >= k; --bit)
            if ((res >> bit) & 1)
                res ^= poly << (bit - k);
        return res & ((1ULL << k) - 1);
    }
}
namespace MathTools
{
    bool isPrime(int n)
    {
        if (n < 2)
            return false;
        for (int i = 2; i * i <= n; i++)
        {
            if (n % i == 0)
                return false;
        }
        return true;
    }
    int getFactor(int n, int start = 2)
    {
        for (int i = start; i * i <= n; i++)
        {
            if (n % i == 0)
                return i;
        }
        return n;
    }
}
namespace OLS
{
    using Square = std::vector<std::vector<int>>;
    using Vi = std::vector<int>;
    std::pair<Square, Square> constructOLS(int n);
    Square makeSquare(int n, int val = 0) { return Square(n, Vi(n, val)); }
    Square formatOA(Square OA)
    {
        if (OA.size() < 2)
            return OA;
        int n = OA[0].size();
        int m = sqrt(n);
        assert(m * m == n);
        Vi p(n);
        for (int j = 0; j < n; j++)
            p[j] = j;
        std::sort(p.begin(), p.end(), [&](int i, int j)
                  { return OA[0][i] ^ OA[0][j] ? OA[0][i] < OA[0][j] : OA[1][i] < OA[1][j]; });
        Square nOA = OA;
        for (int j = 0; j < n; j++)
        {
            for (int i = 0; i < OA.size(); i++)
                OA[i][j] = nOA[i][p[j]];
        }
        Vi vec;
        for (int i = 0; i < OA.size(); i++)
        {
            for (int j = 0; j < OA[i].size(); j++)
                vec.emplace_back(OA[i][j]);
        }
        sort(vec.begin(), vec.end());
        vec.erase(unique(vec.begin(), vec.end()), vec.end());
        for (int i = 0; i < OA.size(); i++)
        {
            for (int j = 0; j < OA[i].size(); j++)
                OA[i][j] = lower_bound(vec.begin(), vec.end(), OA[i][j]) - vec.begin();
        }
        for (int i = 2; i < OA.size(); i++)
        {
            Vi to(m);
            for (int j = 0; j < m; j++)
                to[OA[i][j]] = j;
            for (int j = 0; j < OA[i].size(); j++)
                OA[i][j] = to[OA[i][j]];
        }
        return OA;
    }
    Square genOAbyOLS(std::vector<Square> L)
    {
        if (L.empty())
            return {};
        int len = L[0].size();
        int n = 2 + L.size();
        int m = len * len;
        Square OA(n, Vi(m));
        for (int j = 0; j < m; j++)
            OA[0][j] = j / len;
        for (int j = 0; j < m; j++)
            OA[1][j] = j % len;
        for (int k = 0; k < L.size(); k++)
        {
            for (int i = 0; i < len; i++)
            {
                for (int j = 0; j < len; j++)
                    OA[k + 2][i * len + j] = L[k][i][j];
            }
        }
        OA = formatOA(OA);
        return OA;
    }
    std::pair<Square, Square> genOLSbyOA(Square OA)
    {
        if (OA.size() > 4)
            OA.resize(4);
        if (OA.size() != 4)
            return {{}, {}};
        OA = formatOA(OA);
        int n = OA[0].size();
        int m = sqrt(n);
        assert(m * m == n);
        Square L1 = makeSquare(m, 0);
        Square L2 = makeSquare(m, 0);
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < m; j++)
                L1[i][j] = OA[2][i * m + j];
        }
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < m; j++)
                L2[i][j] = OA[3][i * m + j];
        }
        return {L1, L2};
    }
    std::pair<Square, Square> productOLS(Square A1, Square B1, Square A2, Square B2)
    {
        int n1 = A1.size(), n2 = A2.size();
        int n = n1 * n2;
        Square L1 = makeSquare(n);
        Square L2 = makeSquare(n);
        for (int r = 0; r < n; r++)
        {
            int ra = r / n2, rb = r % n2;
            for (int c = 0; c < n; c++)
            {
                int ca = c / n2, cb = c % n2;
                L1[r][c] = A1[ra][ca] * n2 + A2[rb][cb];
                L2[r][c] = B1[ra][ca] * n2 + B2[rb][cb];
            }
        }
        return {L1, L2};
    }
    std::tuple<Square, Square, Square> productOLS(Square A1, Square B1, Square C1, Square A2, Square B2, Square C2)
    {
        int n1 = A1.size(), n2 = A2.size();
        int n = n1 * n2;
        Square L1 = makeSquare(n);
        Square L2 = makeSquare(n);
        Square L3 = makeSquare(n);
        for (int r = 0; r < n; r++)
        {
            int ra = r / n2, rb = r % n2;
            for (int c = 0; c < n; c++)
            {
                int ca = c / n2, cb = c % n2;
                L1[r][c] = A1[ra][ca] * n2 + A2[rb][cb];
                L2[r][c] = B1[ra][ca] * n2 + B2[rb][cb];
                L3[r][c] = C1[ra][ca] * n2 + C2[rb][cb];
            }
        }
        return {L1, L2, L3};
    }
    std::tuple<Square, Square, Square, Square> constructOLS4(int n)
    {
        if (n == 5)
        {
            Square L1 = makeSquare(n);
            Square L2 = makeSquare(n);
            Square L3 = makeSquare(n);
            Square L4 = makeSquare(n);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    L1[i][j] = (i + 1 * j) % n;
            }
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    L2[i][j] = (i + 2 * j) % n;
            }
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    L3[i][j] = (i + 3 * j) % n;
            }
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    L4[i][j] = (i + 4 * j) % n;
            }
            return {L1, L2, L3, L4};
        }
        throw std::runtime_error("wtf, not implemented");
    }
    std::tuple<Square, Square, Square> constructOLS3(int n)
    {
        if (n == 1)
            return {{{1}}, {{1}}, {{1}}};
        if (n == 2 || n == 6)
            return {{}, {}, {}};
        if (!(n & (n - 1)))
        {
            Square L1 = makeSquare(n);
            Square L2 = makeSquare(n);
            Square L3 = makeSquare(n);
            int k = std::__lg(n);
            uint64_t poly = GF::irreducible_poly(k);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    L1[i][j] = (i ^ GF::gf_mul(1, j, poly, k));
            }
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    L2[i][j] = (i ^ GF::gf_mul(2, j, poly, k));
            }
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    L3[i][j] = (i ^ GF::gf_mul(3, j, poly, k));
            }
            return {L1, L2, L3};
        }
        if (n == 9)
        {
            Square L1 = {{0, 1, 2, 3, 4, 5, 6, 7, 8},
                         {1, 2, 0, 4, 5, 3, 7, 8, 6},
                         {2, 0, 1, 5, 3, 4, 8, 6, 7},
                         {3, 4, 5, 6, 7, 8, 0, 1, 2},
                         {4, 5, 3, 7, 8, 6, 1, 2, 0},
                         {5, 3, 4, 8, 6, 7, 2, 0, 1},
                         {6, 7, 8, 0, 1, 2, 3, 4, 5},
                         {7, 8, 6, 1, 2, 0, 4, 5, 3},
                         {8, 6, 7, 2, 0, 1, 5, 3, 4}};
            Square L2 = {{{0, 1, 2, 3, 4, 5, 6, 7, 8},
                          {6, 7, 8, 0, 1, 2, 3, 4, 5},
                          {3, 4, 5, 6, 7, 8, 0, 1, 2},
                          {1, 2, 0, 4, 5, 3, 7, 8, 6},
                          {7, 8, 6, 1, 2, 0, 4, 5, 3},
                          {4, 5, 3, 7, 8, 6, 1, 2, 0},
                          {2, 0, 1, 5, 3, 4, 8, 6, 7},
                          {8, 6, 7, 2, 0, 1, 5, 3, 4},
                          {5, 3, 4, 8, 6, 7, 2, 0, 1}}};
            Square L3 = {{0, 1, 2, 3, 4, 5, 6, 7, 8},
                         {7, 8, 6, 1, 2, 0, 4, 5, 3},
                         {5, 3, 4, 8, 6, 7, 2, 0, 1},
                         {4, 5, 3, 7, 8, 6, 1, 2, 0},
                         {2, 0, 1, 5, 3, 4, 8, 6, 7},
                         {6, 7, 8, 0, 1, 2, 3, 4, 5},
                         {8, 6, 7, 2, 0, 1, 5, 3, 4},
                         {3, 4, 5, 6, 7, 8, 0, 1, 2},
                         {1, 2, 0, 4, 5, 3, 7, 8, 6}};
            return {L1, L2, L3};
        }
        if (MathTools::isPrime(n))
        // if (n == 5 || n == 13 || n == 19)
        {
            Square L1 = makeSquare(n);
            Square L2 = makeSquare(n);
            Square L3 = makeSquare(n);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    L1[i][j] = (i + 1 * j) % n;
            }
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    L2[i][j] = (i + 2 * j) % n;
            }
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    L3[i][j] = (i + 3 * j) % n;
            }
            return {L1, L2, L3};
        }
        if (n % 4 == 0)
        {
            if (n == 12)
            {
                Square L1 = {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
                             {1, 2, 3, 4, 5, 0, 7, 8, 9, 10, 11, 6},
                             {2, 3, 4, 5, 0, 1, 8, 9, 10, 11, 6, 7},
                             {3, 4, 5, 0, 1, 2, 9, 10, 11, 6, 7, 8},
                             {4, 5, 0, 1, 2, 3, 10, 11, 6, 7, 8, 9},
                             {5, 0, 1, 2, 3, 4, 11, 6, 7, 8, 9, 10},
                             {6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5},
                             {7, 8, 9, 10, 11, 6, 1, 2, 3, 4, 5, 0},
                             {8, 9, 10, 11, 6, 7, 2, 3, 4, 5, 0, 1},
                             {9, 10, 11, 6, 7, 8, 3, 4, 5, 0, 1, 2},
                             {10, 11, 6, 7, 8, 9, 4, 5, 0, 1, 2, 3},
                             {11, 6, 7, 8, 9, 10, 5, 0, 1, 2, 3, 4}};
                Square L2 = {{0, 6, 8, 2, 7, 1, 9, 11, 4, 10, 5, 3},
                             {1, 7, 9, 3, 8, 2, 10, 6, 5, 11, 0, 4},
                             {2, 8, 10, 4, 9, 3, 11, 7, 0, 6, 1, 5},
                             {3, 9, 11, 5, 10, 4, 6, 8, 1, 7, 2, 0},
                             {4, 10, 6, 0, 11, 5, 7, 9, 2, 8, 3, 1},
                             {5, 11, 7, 1, 6, 0, 8, 10, 3, 9, 4, 2},
                             {6, 0, 2, 8, 1, 7, 3, 5, 10, 4, 11, 9},
                             {7, 1, 3, 9, 2, 8, 4, 0, 11, 5, 6, 10},
                             {8, 2, 4, 10, 3, 9, 5, 1, 6, 0, 7, 11},
                             {9, 3, 5, 11, 4, 10, 0, 2, 7, 1, 8, 6},
                             {10, 4, 0, 6, 5, 11, 1, 3, 8, 2, 9, 7},
                             {11, 5, 1, 7, 0, 6, 2, 4, 9, 3, 10, 8}};
                Square L3 = {{0, 3, 6, 1, 9, 11, 2, 8, 5, 4, 7, 10},
                             {1, 4, 7, 2, 10, 6, 3, 9, 0, 5, 8, 11},
                             {2, 5, 8, 3, 11, 7, 4, 10, 1, 0, 9, 6},
                             {3, 0, 9, 4, 6, 8, 5, 11, 2, 1, 10, 7},
                             {4, 1, 10, 5, 7, 9, 0, 6, 3, 2, 11, 8},
                             {5, 2, 11, 0, 8, 10, 1, 7, 4, 3, 6, 9},
                             {6, 9, 0, 7, 3, 5, 8, 2, 11, 10, 1, 4},
                             {7, 10, 1, 8, 4, 0, 9, 3, 6, 11, 2, 5},
                             {8, 11, 2, 9, 5, 1, 10, 4, 7, 6, 3, 0},
                             {9, 6, 3, 10, 0, 2, 11, 5, 8, 7, 4, 1},
                             {10, 7, 4, 11, 1, 3, 6, 0, 9, 8, 5, 2},
                             {11, 8, 5, 6, 2, 4, 7, 1, 10, 9, 0, 3}};
                return {L1, L2, L3};
            }
            else if (n == 24)
            {
                std::vector<Vi> PBD = {{0, 1, 2, 3, 4},
                                       {5, 6, 7, 8, 9},
                                       {10, 11, 12, 13, 14},
                                       {15, 16, 17, 18, 19},
                                       {20, 21, 22, 23},
                                       {0, 5, 10, 15, 20},
                                       {1, 6, 11, 16, 21},
                                       {2, 7, 12, 17, 22},
                                       {3, 8, 13, 18, 23},
                                       {4, 9, 14, 19},
                                       {0, 6, 12, 18},
                                       {1, 7, 13, 19, 20},
                                       {2, 8, 14, 15, 21},
                                       {3, 9, 10, 16, 22},
                                       {4, 5, 11, 17, 23},
                                       {0, 7, 14, 16, 23},
                                       {1, 8, 10, 17},
                                       {2, 9, 11, 18, 20},
                                       {3, 5, 12, 19, 21},
                                       {4, 6, 13, 15, 22},
                                       {0, 8, 11, 19, 22},
                                       {1, 9, 12, 15, 23},
                                       {2, 5, 13, 16},
                                       {3, 6, 14, 17, 20},
                                       {4, 7, 10, 18, 21},
                                       {0, 9, 13, 17, 21},
                                       {1, 5, 14, 18, 22},
                                       {2, 6, 10, 19, 23},
                                       {3, 7, 11, 15},
                                       {4, 8, 12, 16, 20}};
                Square OA(5, Vi(n * n));
                const int q = 4;
                std::map<int, std::vector<Vi>> mp;
                for (int k : {4, 5})
                {
                    if (k == 4)
                    {
                        auto [L1, L2, L3] = constructOLS3(k);
                        auto D = genOAbyOLS({L1, L2, L3});
                        mp[k] = D;
                    }
                    else
                    {
                        auto [L1, L2, L3, L4] = constructOLS4(k);
                        auto D = genOAbyOLS({L1, L2, L3, L4});
                        auto nD = std::vector<Vi>(q + 1, Vi(k * k - k));
                        for (int i = 1; i < q + 2; i++)
                        {
                            for (int j = k; j < k * k; j++)
                                nD[i - 1][j - k] = D[i][j];
                        }
                        mp[k] = nD;
                    }
                }
                int cnt = 0;
                std::set<int> occ;
                for (auto B : PBD)
                {
                    int k = B.size();
                    for (int x : B)
                        occ.insert(x);
                    assert(mp.count(k));
                    auto P = mp[k];
                    std::map<int, int> f;
                    int id = 0;
                    for (int i = 0; i < P.size(); i++)
                    {
                        for (int j = 0; j < P[i].size(); j++)
                        {
                            if (f.count(P[i][j]))
                                P[i][j] = f[P[i][j]];
                            else
                                P[i][j] = (f[P[i][j]] = B[id++]);
                        }
                    }
                    int len = P[0].size();
                    for (int j = 0; j < len; j++)
                    {
                        for (int i = 0; i < q + 1; i++)
                            OA[i][cnt] = P[i][j];
                        cnt++;
                    }
                }
                for (auto B : PBD)
                {
                    int k = B.size();
                    if (k == 4)
                    {
                        for (int x : B)
                            occ.erase(x);
                    }
                }
                Square F(5, Vi(occ.size()));
                int tc = 0;
                for (int j : occ)
                {
                    for (int i = 0; i < q + 1; i++)
                        F[i][tc] = j;
                    tc++;
                }
                int len = F[0].size();
                for (int j = 0; j < len; j++)
                {
                    for (int i = 0; i < q + 1; i++)
                        OA[i][cnt] = F[i][j];
                    cnt++;
                }
                OA = formatOA(OA);
                int n = OA[0].size();
                int m = sqrt(n);
                assert(m * m == n);
                Square L1 = makeSquare(m, 0);
                Square L2 = makeSquare(m, 0);
                Square L3 = makeSquare(m, 0);
                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < m; j++)
                        L1[i][j] = OA[2][i * m + j];
                }
                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < m; j++)
                        L2[i][j] = OA[3][i * m + j];
                }
                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < m; j++)
                        L3[i][j] = OA[4][i * m + j];
                }
                return {L1, L2, L3};
            }
            else
            {
                if (n % 9 == 0)
                {
                    auto [A1, B1, C1] = constructOLS3(9);
                    auto [A2, B2, C2] = constructOLS3(n / 9);
                    return productOLS(A1, B1, C1, A2, B2, C2);
                }
                if (n % 3 == 0)
                {
                    int pw2 = n & (-n);
                    int lg = __builtin_ctz(pw2);
                    if (lg & 1)
                    {
                        auto [A1, B1, C1] = constructOLS3(24);
                        auto [A2, B2, C2] = constructOLS3(n / 24);
                        return productOLS(A1, B1, C1, A2, B2, C2);
                    }
                    else
                    {
                        auto [A1, B1, C1] = constructOLS3(12);
                        auto [A2, B2, C2] = constructOLS3(n / 12);
                        return productOLS(A1, B1, C1, A2, B2, C2);
                    }
                }
                int pw2 = n & (-n);
                int lg = __builtin_ctz(pw2);
                if (lg & 1)
                {
                    auto [A1, B1, C1] = constructOLS3(8);
                    auto [A2, B2, C2] = constructOLS3(n / 8);
                    return productOLS(A1, B1, C1, A2, B2, C2);
                }
                else
                {
                    auto [A1, B1, C1] = constructOLS3(4);
                    auto [A2, B2, C2] = constructOLS3(n / 4);
                    return productOLS(A1, B1, C1, A2, B2, C2);
                }
            }
        }
        int d = MathTools::getFactor(n);
        auto [A1, B1, C1] = constructOLS3(d);
        auto [A2, B2, C2] = constructOLS3(n / d);
        return productOLS(A1, B1, C1, A2, B2, C2);
        throw std::runtime_error("wtf, not implemented");
        return {{}, {}, {}};
    }
    std::pair<Square, Square> constructOLSbyPBD(int n, std::vector<Vi> PBD, std::set<int> K0, std::set<int> K)
    {
        Square OA(4, Vi(n * n));
        const int q = 3;
        std::map<int, std::vector<Vi>> mp;
        for (int k : K)
        {
            if (K0.count(k))
            {
                auto [L1, L2] = constructOLS(k);
                auto D = genOAbyOLS({L1, L2});
                mp[k] = D;
            }
            else
            {
                auto [L1, L2, L3] = constructOLS3(k);
                auto D = genOAbyOLS({L1, L2, L3});
                auto nD = std::vector<Vi>(q + 1, Vi(k * k - k));
                for (int i = 1; i < q + 2; i++)
                {
                    for (int j = k; j < k * k; j++)
                        nD[i - 1][j - k] = D[i][j];
                }
                mp[k] = nD;
            }
        }
        int cnt = 0;
        std::set<int> occ;
        for (auto B : PBD)
        {
            int k = B.size();
            for (int x : B)
                occ.insert(x);
            assert(mp.count(k));
            auto P = mp[k];
            std::map<int, int> f;
            int id = 0;
            for (int i = 0; i < P.size(); i++)
            {
                for (int j = 0; j < P[i].size(); j++)
                {
                    if (f.count(P[i][j]))
                        P[i][j] = f[P[i][j]];
                    else
                        P[i][j] = (f[P[i][j]] = B[id++]);
                }
            }
            int len = P[0].size();
            for (int j = 0; j < len; j++)
            {
                for (int i = 0; i < q + 1; i++)
                    OA[i][cnt] = P[i][j];
                cnt++;
            }
        }
        for (auto B : PBD)
        {
            int k = B.size();
            if (K0.count(k))
            {
                for (int x : B)
                    occ.erase(x);
            }
        }
        Square F(4, Vi(occ.size()));
        int tc = 0;
        for (int j : occ)
        {
            for (int i = 0; i < q + 1; i++)
                F[i][tc] = j;
            tc++;
        }
        int len = F[0].size();
        for (int j = 0; j < len; j++)
        {
            for (int i = 0; i < q + 1; i++)
                OA[i][cnt] = F[i][j];
            cnt++;
        }
        return genOLSbyOA(OA);
    }
    std::pair<std::vector<std::vector<Vi>>, std::vector<std::vector<Vi>>> getParallelPBD(int m, int k) // k = 4
    {
        assert(k == 4);
        auto [L1, L2, L3] = constructOLS3(m);
        auto OA = genOAbyOLS({L1, L2, L3});
        OA.erase(OA.begin());
        int t = 0;
        for (int i = 0; i < OA.size(); i++)
        {
            for (int j = 0; j < OA[i].size(); j++)
                OA[i][j] += t;
            t += m;
        }
        std::vector<std::vector<Vi>> PBD_L;
        std::vector<std::vector<Vi>> PBD_R;
        for (int l = 0, r = m; l < m * m; l += m, r += m)
        {
            std::vector<Vi> vec;
            for (int j = l; j < r; j++)
            {
                Vi tg;
                for (int i = 0; i < k; i++)
                    tg.emplace_back(OA[i][j]);
                vec.emplace_back(tg);
            }
            PBD_L.emplace_back(vec);
        }
        int l = 0, r = m;
        {
            std::vector<Vi> vec;
            for (int i = 0; i < k; i++)
            {
                Vi tg;
                for (int j = l; j < r; j++)
                    tg.emplace_back(OA[i][j]);
                vec.emplace_back(tg);
            }
            PBD_R.emplace_back(vec);
        }
        return {PBD_L, PBD_R};
    }
    std::pair<Square, Square> constructMKX(int m, int k, int x)
    {
        auto [PBD_L, PBD_R] = getParallelPBD(m, k);
        std::vector<Vi> PBD;
        int cnt = m * k;
        for (auto B : PBD_L)
        {
            for (auto v : B)
            {
                if (cnt < m * k + x)
                    v.emplace_back(cnt);
                PBD.emplace_back(v);
            }
            cnt++;
        }
        for (auto B : PBD_R)
        {
            for (auto v : B)
                PBD.emplace_back(v);
        }
        Vi Bs;
        for (int i = m * k; i < m * k + x; i++)
            Bs.emplace_back(i);
        PBD.emplace_back(Bs);
        int n = m * k + x;
        return constructOLSbyPBD(n, PBD, {x}, {x, m, k, k + 1});
    }
    std::pair<Square, Square> constructOLS(int n)
    {
        if (n == 2 || n == 6)
            throw std::runtime_error("No pair of OLS exists for n = 2 or n = 6");
        if (n % 2 == 1) // 2ez bro
        {
            Square L1 = makeSquare(n, 0);
            Vi vec(n);
            for (int i = 0; i < n; i++)
                vec[i] = i;
            for (int i = 0; i < n; i++)
            {
                L1[i] = vec;
                std::rotate(vec.begin(), next(vec.begin()), vec.end());
            }
            Square L2 = L1;
            for (int i = 0; i < n; i++)
                reverse(L2[i].begin(), L2[i].end());
            return {L1, L2};
        }
        if (!(n & (n - 1)))
        {
            auto [L1, L2, L3] = constructOLS3(n);
            return {L1, L2};
        }
        if (n % 4 == 0)
        {
            for (int d = 4; d <= 8; d += 4)
            {
                if (n % d == 0 && n / d != 2 && n / d != 6)
                {
                    auto [A1, B1] = constructOLS(d);
                    auto [A2, B2] = constructOLS(n / d);
                    return productOLS(A1, B1, A2, B2);
                }
            }
            return {{}, {}};
        }
        // very hard: n = 2 (mod 4)
        if (n == 10)
        {
            /*
            thanks to R. C. BOSE, S. S. SHRIKHANDE, AND E. T. PARKER
            https://www.cambridge.org/core/services/aop-cambridge-core/content/view/1152262BF046F8632638FE9C10610136/S0008414X0000986Xa.pdf/further-results-on-the-construction-of-mutually-orthogonal-latin-squares-and-the-falsity-of-eulers-conjecture.pdf
            */
            Square L1 = {{0, 6, 5, 4, 9, 8, 7, 1, 2, 3},
                         {7, 1, 0, 6, 5, 9, 8, 2, 3, 4},
                         {8, 7, 2, 1, 0, 6, 9, 3, 4, 5},
                         {9, 8, 7, 3, 2, 1, 0, 4, 5, 6},
                         {1, 9, 8, 7, 4, 3, 2, 5, 6, 0},
                         {3, 2, 9, 8, 7, 5, 4, 6, 0, 1},
                         {5, 4, 3, 9, 8, 7, 6, 0, 1, 2},
                         {2, 3, 4, 5, 6, 0, 1, 7, 8, 9},
                         {4, 5, 6, 0, 1, 2, 3, 8, 9, 7},
                         {6, 0, 1, 2, 3, 4, 5, 9, 7, 8}};
            Square L2 = {{0, 7, 8, 9, 1, 3, 5, 2, 4, 6},
                         {6, 1, 7, 8, 9, 2, 4, 3, 5, 0},
                         {5, 0, 2, 7, 8, 9, 3, 4, 6, 1},
                         {4, 6, 1, 3, 7, 8, 9, 5, 0, 2},
                         {9, 5, 0, 2, 4, 7, 8, 6, 1, 3},
                         {8, 9, 6, 1, 3, 5, 7, 0, 2, 4},
                         {7, 8, 9, 0, 2, 4, 6, 1, 3, 5},
                         {1, 2, 3, 4, 5, 6, 0, 7, 8, 9},
                         {2, 3, 4, 5, 6, 0, 1, 9, 7, 8},
                         {3, 4, 5, 6, 0, 1, 2, 8, 9, 7}};
            return {L1, L2};
        }
        if (n == 18)
        {
            return constructOLSbyPBD(18,
                                     {{1, 4, 14, 16},
                                      {1, 15, 17},
                                      {3, 6, 16, 18},
                                      {3, 4, 7, 17, 19},
                                      {4, 8, 18, 20},
                                      {6, 9, 19},
                                      {1, 6, 7, 10, 20},
                                      {7, 8, 11},
                                      {1, 3, 8, 9, 12},
                                      {4, 9, 10, 13},
                                      {3, 10, 11, 14},
                                      {4, 6, 11, 12, 15},
                                      {7, 12, 13, 16},
                                      {6, 8, 13, 14, 17},
                                      {7, 9, 14, 15, 18},
                                      {8, 10, 15, 16, 19},
                                      {9, 11, 16, 17, 20},
                                      {10, 12, 17, 18},
                                      {1, 11, 13, 18, 19},
                                      {12, 14, 19, 20},
                                      {3, 13, 15, 20}},
                                     {3}, {3, 4, 5});
        }
        if (n == 14)
        {
            /*
            thanks to MatrixGroup.
            https://www.luogu.com.cn/article/i8x6is6u
            */
            Square L1 = {{0, 8, 3, 12, 9, 2, 5, 10, 6, 11, 1, 4, 13, 7},
                         {13, 1, 9, 4, 0, 10, 3, 6, 11, 7, 12, 2, 5, 8},
                         {6, 13, 2, 10, 5, 1, 11, 4, 7, 12, 8, 0, 3, 9},
                         {4, 7, 13, 3, 11, 6, 2, 12, 5, 8, 0, 9, 1, 10},
                         {2, 5, 8, 13, 4, 12, 7, 3, 0, 6, 9, 1, 10, 11},
                         {11, 3, 6, 9, 13, 5, 0, 8, 4, 1, 7, 10, 2, 12},
                         {3, 12, 4, 7, 10, 13, 6, 1, 9, 5, 2, 8, 11, 0},
                         {12, 4, 0, 5, 8, 11, 13, 7, 2, 10, 6, 3, 9, 1},
                         {10, 0, 5, 1, 6, 9, 12, 13, 8, 3, 11, 7, 4, 2},
                         {5, 11, 1, 6, 2, 7, 10, 0, 13, 9, 4, 12, 8, 3},
                         {9, 6, 12, 2, 7, 3, 8, 11, 1, 13, 10, 5, 0, 4},
                         {1, 10, 7, 0, 3, 8, 4, 9, 12, 2, 13, 11, 6, 5},
                         {7, 2, 11, 8, 1, 4, 9, 5, 10, 0, 3, 13, 12, 6},
                         {8, 9, 10, 11, 12, 0, 1, 2, 3, 4, 5, 6, 7, 13}};
            Square L2 = makeSquare(n, 0);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    L2[i][j] = L1[j][i];
            }
            return {L1, L2};
        }
        if (n == 22)
        {
            /*
            thanks to BY R. C. BOSE AND S. S. SHRIKHANDE
            https://www.pnas.org/doi/epdf/10.1073/pnas.45.5.734
            */
            Square L1 = {{0, 3, 6, 15, 5, 19, 21, 14, 18, 20, 11, 17, 9, 8, 16, 1, 7, 10, 13, 4, 12, 2},
                         {15, 1, 4, 0, 16, 6, 20, 9, 14, 19, 21, 12, 18, 10, 17, 3, 2, 8, 11, 7, 5, 13},
                         {21, 16, 2, 5, 1, 17, 0, 11, 10, 14, 20, 15, 13, 19, 18, 7, 4, 3, 9, 12, 8, 6},
                         {1, 15, 17, 3, 6, 2, 18, 20, 12, 11, 14, 21, 16, 7, 19, 0, 8, 5, 4, 10, 13, 9},
                         {19, 2, 16, 18, 4, 0, 3, 8, 21, 13, 12, 14, 15, 17, 20, 10, 1, 9, 6, 5, 11, 7},
                         {4, 20, 3, 17, 19, 5, 1, 18, 9, 15, 7, 13, 14, 16, 21, 8, 11, 2, 10, 0, 6, 12},
                         {2, 5, 21, 4, 18, 20, 6, 17, 19, 10, 16, 8, 7, 14, 15, 13, 9, 12, 3, 11, 1, 0},
                         {16, 19, 15, 13, 21, 10, 12, 7, 4, 1, 18, 2, 17, 20, 0, 11, 14, 6, 5, 9, 3, 8},
                         {13, 17, 20, 16, 7, 15, 11, 21, 8, 5, 2, 19, 3, 18, 1, 9, 12, 14, 0, 6, 10, 4},
                         {12, 7, 18, 21, 17, 8, 16, 19, 15, 9, 6, 3, 20, 4, 2, 5, 10, 13, 14, 1, 0, 11},
                         {17, 13, 8, 19, 15, 18, 9, 5, 20, 16, 10, 0, 4, 21, 3, 12, 6, 11, 7, 14, 2, 1},
                         {10, 18, 7, 9, 20, 16, 19, 15, 6, 21, 17, 11, 1, 5, 4, 2, 13, 0, 12, 8, 14, 3},
                         {20, 11, 19, 8, 10, 21, 17, 6, 16, 0, 15, 18, 12, 2, 5, 4, 3, 7, 1, 13, 9, 14},
                         {18, 21, 12, 20, 9, 11, 15, 3, 0, 17, 1, 16, 19, 13, 6, 14, 5, 4, 8, 2, 7, 10},
                         {7, 8, 9, 10, 11, 12, 13, 16, 17, 18, 19, 20, 21, 15, 14, 6, 0, 1, 2, 3, 4, 5},
                         {3, 0, 11, 1, 12, 9, 14, 2, 5, 8, 4, 7, 10, 6, 13, 15, 17, 19, 21, 16, 18, 20},
                         {14, 4, 1, 12, 2, 13, 10, 0, 3, 6, 9, 5, 8, 11, 7, 21, 16, 18, 20, 15, 17, 19},
                         {11, 14, 5, 2, 13, 3, 7, 12, 1, 4, 0, 10, 6, 9, 8, 20, 15, 17, 19, 21, 16, 18},
                         {8, 12, 14, 6, 3, 7, 4, 10, 13, 2, 5, 1, 11, 0, 9, 19, 21, 16, 18, 20, 15, 17},
                         {5, 9, 13, 14, 0, 4, 8, 1, 11, 7, 3, 6, 2, 12, 10, 18, 20, 15, 17, 19, 21, 16},
                         {9, 6, 10, 7, 14, 1, 5, 13, 2, 12, 8, 4, 0, 3, 11, 17, 19, 21, 16, 18, 20, 15},
                         {6, 10, 0, 11, 8, 14, 2, 4, 7, 3, 13, 9, 5, 1, 12, 16, 18, 20, 15, 17, 19, 21}};
            Square L2 = {{0, 15, 21, 1, 19, 4, 2, 16, 13, 12, 17, 10, 20, 18, 7, 3, 14, 11, 8, 5, 9, 6},
                         {3, 1, 16, 15, 2, 20, 5, 19, 17, 7, 13, 18, 11, 21, 8, 0, 4, 14, 12, 9, 6, 10},
                         {6, 4, 2, 17, 16, 3, 21, 15, 20, 18, 8, 7, 19, 12, 9, 11, 1, 5, 14, 13, 10, 0},
                         {15, 0, 5, 3, 18, 17, 4, 13, 16, 21, 19, 9, 8, 20, 10, 1, 12, 2, 6, 14, 7, 11},
                         {5, 16, 1, 6, 4, 19, 18, 21, 7, 17, 15, 20, 10, 9, 11, 12, 2, 13, 3, 0, 14, 8},
                         {19, 6, 17, 2, 0, 5, 20, 10, 15, 8, 18, 16, 21, 11, 12, 9, 13, 3, 7, 4, 1, 14},
                         {21, 20, 0, 18, 3, 1, 6, 12, 11, 16, 9, 19, 17, 15, 13, 14, 10, 7, 4, 8, 5, 2},
                         {14, 9, 11, 20, 8, 18, 17, 7, 21, 19, 5, 15, 6, 3, 16, 2, 0, 12, 10, 1, 13, 4},
                         {18, 14, 10, 12, 21, 9, 19, 4, 8, 15, 20, 6, 16, 0, 17, 5, 3, 1, 13, 11, 2, 7},
                         {20, 19, 14, 11, 13, 15, 10, 1, 5, 9, 16, 21, 0, 17, 18, 8, 6, 4, 2, 7, 12, 3},
                         {11, 21, 20, 14, 12, 7, 16, 18, 2, 6, 10, 17, 15, 1, 19, 4, 9, 0, 5, 3, 8, 13},
                         {17, 12, 15, 21, 14, 13, 8, 2, 19, 3, 0, 11, 18, 16, 20, 7, 5, 10, 1, 6, 4, 9},
                         {9, 18, 13, 16, 15, 14, 7, 17, 3, 20, 4, 1, 12, 19, 21, 10, 8, 6, 11, 2, 0, 5},
                         {8, 10, 19, 7, 17, 16, 14, 20, 18, 4, 21, 5, 2, 13, 15, 6, 11, 9, 0, 12, 3, 1},
                         {16, 17, 18, 19, 20, 21, 15, 0, 1, 2, 3, 4, 5, 6, 14, 13, 7, 8, 9, 10, 11, 12},
                         {1, 3, 7, 0, 10, 8, 13, 11, 9, 5, 12, 2, 4, 14, 6, 15, 19, 16, 20, 17, 21, 18},
                         {7, 2, 4, 8, 1, 11, 9, 14, 12, 10, 6, 13, 3, 5, 0, 19, 16, 20, 17, 21, 18, 15},
                         {10, 8, 3, 5, 9, 2, 12, 6, 14, 13, 11, 0, 7, 4, 1, 16, 20, 17, 21, 18, 15, 19},
                         {13, 11, 9, 4, 6, 10, 3, 5, 0, 14, 7, 12, 1, 8, 2, 20, 17, 21, 18, 15, 19, 16},
                         {4, 7, 12, 10, 5, 0, 11, 9, 6, 1, 14, 8, 13, 2, 3, 17, 21, 18, 15, 19, 16, 20},
                         {12, 5, 8, 13, 11, 6, 1, 3, 10, 0, 2, 14, 9, 7, 4, 21, 18, 15, 19, 16, 20, 17},
                         {2, 13, 6, 9, 7, 12, 0, 8, 4, 11, 1, 3, 14, 10, 5, 18, 15, 19, 16, 20, 17, 21}};
            return {L1, L2};
        }
        if (n == 26)
        {
            Vi W = {0, 19, 20, 21, 22, 23, 24, 25};
            Vi X = {3, 15, 10, 7, 8, 12, 9, 6};
            Vi Y = {1, 0, 0, 0, 0, 0, 0, 0};
            Vi Z = {6, 1, 2, 4, 6, 7, 8, 10};
            Square A;
            A.resize(4);
            auto bindVi = [&](std::vector<Vi> A)
            {
                Vi ans;
                for (auto V : A)
                {
                    for (auto v : V)
                        ans.emplace_back(v);
                }
                return ans;
            };
            A[0] = bindVi({W, Z, X, Y});
            A[1] = bindVi({X, Y, W, Z});
            A[2] = bindVi({Y, W, Z, X});
            A[3] = bindVi({Z, X, Y, W});
            Square OA(4, Vi(n * n));
            int cnt = 0;
            for (int j = 0; j < 19; j++)
            {
                for (int i = 0; i < 4; i++)
                    OA[i][cnt] = j;
                cnt++;
            }
            for (int k = 0; k < 19; k++)
            {
                for (int j = 0; j < 8 * 4; j++)
                {
                    for (int i = 0; i < 4; i++)
                    {
                        if (A[i][j] > 18)
                            OA[i][cnt] = A[i][j];
                        else
                            OA[i][cnt] = (A[i][j] + k) % 19;
                    }
                    cnt++;
                }
            }
            auto [L1, L2] = constructOLS(7);
            Square As = genOAbyOLS({L1, L2});
            for (int j = 0; j < 7 * 7; j++)
            {
                for (int i = 0; i < 4; i++)
                    OA[i][cnt] = As[i][j] + 19;
                cnt++;
            }
            return genOLSbyOA(OA);
        }
        if (n == 38)
        {
            std::set<std::set<int>> BIBD;
            for (int i = 0; i < 41; i++)
            {
                BIBD.insert({(0 + i) % 41, (7 + i) % 41, (10 + i) % 41, (11 + i) % 41, (23 + i) % 41});
                BIBD.insert({(0 + i) % 41, (5 + i) % 41, (14 + i) % 41, (20 + i) % 41, (22 + i) % 41});
            }
            int a = 0, b = 1, c = 2; // found by brute force
            Square PBD;
            for (int i = 0; i < 41; i++)
            {
                std::set<int> s = {(0 + i) % 41, (7 + i) % 41, (10 + i) % 41, (11 + i) % 41, (23 + i) % 41};
                s.erase(a), s.erase(b), s.erase(c);
                PBD.emplace_back(Vi(s.begin(), s.end()));
                s = {(0 + i) % 41, (5 + i) % 41, (14 + i) % 41, (20 + i) % 41, (22 + i) % 41};
                s.erase(a), s.erase(b), s.erase(c);
                PBD.emplace_back(Vi(s.begin(), s.end()));
            }
            return constructOLSbyPBD(n, PBD, {3}, {3, 4, 5});
        }
        if ((n - 1) % 3 == 0) // n = 3m + 1
        {
            int m = (n - 1) / 3;
            Vi W;
            for (int i = 0; i < m; i++)
                W.emplace_back(0);
            Vi X;
            for (int i = 0; i < m; i++)
                X.emplace_back(i + 1);
            Vi Y;
            for (int i = 0; i < m; i++)
                Y.emplace_back(2 * m - i);
            Vi Z;
            for (int i = 0; i < m; i++)
                Z.emplace_back(-(i + 1));
            Square A;
            A.resize(4);
            auto bindVi = [&](std::vector<Vi> A)
            {
                Vi ans;
                for (auto V : A)
                {
                    for (auto v : V)
                        ans.emplace_back(v);
                }
                return ans;
            };
            A[0] = bindVi({W, X, Y, Z});
            A[1] = bindVi({X, W, Z, Y});
            A[2] = bindVi({Y, Z, W, X});
            A[3] = bindVi({Z, Y, X, W});
            Square OA(4, Vi(n * n));
            int cnt = 0;
            for (int j = 0; j < 2 * m + 1; j++)
            {
                for (int i = 0; i < 4; i++)
                    OA[i][cnt] = j;
                cnt++;
            }
            for (int k = 0; k < 2 * m + 1; k++)
            {
                for (int j = 0; j < m * 4; j++)
                {
                    for (int i = 0; i < 4; i++)
                    {
                        if (A[i][j] < 0)
                            OA[i][cnt] = 2 * m - A[i][j];
                        else
                            OA[i][cnt] = (A[i][j] + k) % (2 * m + 1);
                    }
                    cnt++;
                }
            }
            auto [L1, L2] = constructOLS(m);
            Square As = genOAbyOLS({L1, L2});
            for (int j = 0; j < m * m; j++)
            {
                for (int i = 0; i < 4; i++)
                    OA[i][cnt] = As[i][j] + (2 * m + 1);
                cnt++;
            }
            return genOLSbyOA(OA);
        }
        if (n == 62 || n == 86)
        {
            int m = n == 62 ? 13 : 19, k = 4, x = 10;
            return constructMKX(m, k, x);
        }
        for (int d = 2; d < n; d++)
        {
            if (d == 2 || d == 6)
                continue;
            if (n % d == 0 && n / d != 2 && n / d != 6)
            {
                auto [A1, B1] = constructOLS(d);
                auto [A2, B2] = constructOLS(n / d);
                return productOLS(A1, B1, A2, B2);
            }
        }
        if (n % 4 == 2)
        {
            int x = n % 16 == 2 ? 18 : n % 16 == 6 ? 22
                                   : n % 16 == 10  ? 10
                                                   : 14;
            int m = (n - x) / 4, k = 4;
            return constructMKX(m, k, x);
        }
        return {{}, {}};
    }
    bool checkPermutation(Vi p, int n)
    {
        if (p.size() != n)
            return false;
        Vi vis(n);
        for (int i = 0; i < n; i++)
        {
            if (!(0 <= p[i] && p[i] < n))
                return false;
            if (vis[p[i]])
                return false;
            vis[p[i]] = true;
        }
        return true;
    }
    bool checkLS(Square L)
    {
        int n = L.size();
        for (int i = 0; i < n; i++)
        {
            if (!checkPermutation(L[i], n))
                return false;
        }
        for (int j = 0; j < n; j++)
        {
            Vi tmp;
            for (int i = 0; i < n; i++)
                tmp.emplace_back(L[i][j]);
            if (!checkPermutation(tmp, n))
                return false;
        }
        return true;
    }
    bool checkOLS(Square L1, Square L2, int n)
    {
        if (L1.size() != n || L2.size() != n)
            return false;
        if (!checkLS(L1) || !checkLS(L2))
            return false;
        Square vis = makeSquare(n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (vis[L1[i][j]][L2[i][j]])
                    return false;
                vis[L1[i][j]][L2[i][j]] = 1;
            }
        }
        return true;
    }
}; // namespace OLS
#endif // PROJECT_OLS_H_