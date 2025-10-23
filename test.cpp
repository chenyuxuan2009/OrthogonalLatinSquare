#include "ols.h"
#include <iostream>
#include <vector>
int main()
{
    if (0)
    {
        // std::vector<std::vector<int>> OA = {{1, 1, 1, 2, 2, 2, 3, 3, 3},
        //                                     {1, 2, 3, 1, 2, 3, 1, 2, 3},
        //                                     {1, 2, 3, 3, 1, 2, 2, 3, 1},
        //                                     {1, 2, 3, 2, 3, 1, 3, 1, 2}};
        std::vector<std::vector<int>> OA = {{1, 2, 3, 2, 3, 1, 3, 1, 2},
                                            {1, 2, 3, 3, 1, 2, 2, 3, 1},
                                            {1, 1, 1, 2, 2, 2, 3, 3, 3},
                                            {1, 2, 3, 1, 2, 3, 1, 2, 3}};
        auto nOA = OLS::formatOA(OA);
        for (int i = 0; i < nOA.size(); i++)
        {
            for (int j = 0; j < nOA[i].size(); j++)
                printf("%d ", nOA[i][j]);
            printf("\n");
        }
        printf("\n");
        auto [L1, L2] = OLS::genOLSbyOA(OA);
        for (int i = 0; i < L1.size(); i++)
        {
            for (int j = 0; j < L1[i].size(); j++)
                printf("%d ", L1[i][j]);
            printf("\n");
        }
        printf("\n");
        for (int i = 0; i < L2.size(); i++)
        {
            for (int j = 0; j < L2[i].size(); j++)
                printf("%d ", L2[i][j]);
            printf("\n");
        }
        printf("\n");
    }
    if (0)
    {
        auto OLS_3_check = [&](int n)
        {
            auto [L1, L2, L3] = OLS::constructOLS3(n);
            assert(OLS::checkOLS(L1, L2, n));
            assert(OLS::checkOLS(L1, L3, n));
            assert(OLS::checkOLS(L2, L3, n));
            for (int i = 0; i < L1.size(); i++)
            {
                for (int j = 0; j < L1[i].size(); j++)
                    printf("%d ", L1[i][j]);
                printf("\n");
            }
            printf("\n");
            for (int i = 0; i < L2.size(); i++)
            {
                for (int j = 0; j < L2[i].size(); j++)
                    printf("%d ", L2[i][j]);
                printf("\n");
            }
            printf("\n");
            for (int i = 0; i < L3.size(); i++)
            {
                for (int j = 0; j < L3[i].size(); j++)
                    printf("%d ", L3[i][j]);
                printf("\n");
            }
            printf("\n");
        };
        OLS_3_check(4);
        OLS_3_check(5);
        OLS_3_check(8);
        OLS_3_check(9);
        OLS_3_check(12);
        OLS_3_check(24);
    }
    if (1)
    {
        std::vector<int> failTests;
        for (int n = 1; n <= 500; n++)
        {
            if (n == 2 || n == 6) // n = 2 and n = 6 are not orthogonal
                continue;
            auto [L1, L2] = OLS::constructOLS(n);
            if (OLS::checkOLS(L1, L2, n))
                printf("n = %d: orthogonal passed\n", n);
            else
                printf("n = %d: orthogonal failed\n", n), failTests.emplace_back(n);
        }
        if (failTests.empty())
            printf("All tests passed\n");
        else
        {
            printf("Failed tests: ");
            for (int n : failTests)
                printf("%d ", n);
            printf("\n");
        }
    }
    return 0;
}