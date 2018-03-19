#include "gtest/gtest.h"
#include "PCMC/Vector.h"
#include "stdio.h"

namespace PCA{

TEST(Vector, Constructor){
    Vector vector(1.1, 2.2, 3.3);
    EXPECT_DOUBLE_EQ(vector.x, 1.1) << "init vector.x faild";
    EXPECT_DOUBLE_EQ(vector.y, 2.2);
    EXPECT_DOUBLE_EQ(vector.z, 3.3);
}

TEST(Vector, Eqauality){
    Vector vector1(1.1, 2.2, 3.3);
    Vector vector2(1.1, 2.2, 3.3);
    
    char* err;
    err = new char[100];
    sprintf(err, "%.16le", PCA_NUMERICAL_ERROR);
    
    
    EXPECT_TRUE(vector1 == vector2) << err;
    delete [] err;
}

TEST(Macros, equality){
    double a = 10.01010101;
    double b = 10.01010101;
    EXPECT_TRUE(_PCA_IS_EQUAL(a, b));
    
    a = 10.010101010;
    b = 10.010101011;
    EXPECT_FALSE(_PCA_IS_EQUAL(a, b));
}

TEST(Macros, voidPointer){
    double* ptr = nullptr;
    
    EXPECT_EXIT(_PCA_CATCH_VOID_POINTER(ptr, "test"), ::testing::ExitedWithCode(1), "");
}

}//namespace
