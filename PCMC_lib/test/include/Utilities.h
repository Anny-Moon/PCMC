#include "gtest/gtest.h"
#include "PCMC/Utilities.h"

namespace PCA{

TEST(Utilities, MeanValue){
    
    double array [] = {1.1, 2.2, 3.3};
    EXPECT_DOUBLE_EQ(meanValue(3, array), 2.2);

    std::vector<double> vector = {1.1, 2.2, 3.3};
    EXPECT_DOUBLE_EQ(meanValue(3, array), 2.2);
}

TEST(Utilities, commonDivisor){
    EXPECT_EQ(commonDivisor(10,10), 10);
    EXPECT_EQ(commonDivisor(24,12), 12);
    EXPECT_EQ(commonDivisor(12,24), 12);
    EXPECT_EQ(commonDivisor(64,24), 8);
    EXPECT_EQ(commonDivisor(64,24,7), 4); // less or equal than 7
    EXPECT_EQ(commonDivisor(7,9), 1);
}

//TEST(Macros, voidPointer){
//    double* ptr = nullptr;
    
//    EXPECT_EXIT(_PCA_CATCH_VOID_POINTER(ptr, "test"), ::testing::ExitedWithCode(1), "");
//}


}//namespace
