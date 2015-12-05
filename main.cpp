#include <iostream>
#include <stdio.h>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <queue>
#include <fstream>
#include <sstream>
#define random(x) (rand()%x)
#define ROW_ONES_NUM 10000
using namespace std;
typedef unsigned char  U8,  *PU8;
typedef unsigned short U16, *PU16;
typedef unsigned int   U32, *PU32;
using std::vector;
using std::queue;
class CData
{
private:
  U8* m_Data;
  U32  m_Len;
public:

  CData(void): m_Data(NULL), m_Len(0)
  { 
  }
  ~CData(void)
  {    
  }
  CData(U8* data, U32 len): m_Data(NULL), m_Len(0) {SetData(data, len);}; 
  U8 *GetData(void) const { return m_Data; }
  U32 GetLen(void) const { return m_Len; }
  void FreeData()
  {
    if (m_Data)
    {
      delete[] m_Data;
      m_Data = NULL;
    }
  }

  void SetData(CData* data)
  {
    return SetData(data->GetData(), data->GetLen());
  }
  void SetData(const U8* data, U32 len)
  {
    FreeData();
    if (len > 0 && data)
    {
      m_Data = new U8[len];
      memcpy(m_Data, data, len);
      m_Len = len;
    }
  }
  void XorData(CData* data)
  {
    return XorData(data->GetData(), data->GetLen());
  }

  void XorData(const U8* data, U32 len)
  {
    if (len > m_Len)
    {
      U8* tmp = new U8[len];
      memset(tmp,0,len);
      if (m_Len > 0)
      {
	memcpy(tmp, m_Data, m_Len);
      }
      FreeData();
      m_Data = tmp;
      m_Len = len;
    }
    
    for (U32 i = 0; i < m_Len; ++i)
    {
      m_Data[i] ^= data[i];
    }
  }

};
void decode(U32 M, U32 L, U8** A, CData** D,CData** C,U16* m_Degree)
{

  U32* c = new U32[L];
  U32* d = new U32[M];
  for (U32 i = 0; i < M; ++i)
  {
    if (i < L)
    {
      c[i] = i;
    }
    d[i] = i;
  }
  for (U32 i = 0; i < M; ++i)
  {
    U32 length = (*D[d[i]]).GetLen();
    for (U32 j = 0; j < length; ++j)
    {
      printf("%6d", (*D[d[i]]).GetData()[j]);
    }
    if (length==0)
      printf(" 0 ");
  }
  /*-----------------------------------First Phase------------------------------------*/
  U32 zerosRowNum = 0;
  U32 i = 0; // i is the rank of Identity Matrix I
  U32 u = 0; // u is the column number of Matrix U
  while (i + u < L)
  {
    //Let r be the minimum integer such that at least one row of A has
    //exactly r ones in V.
    U32 ones_pos[ROW_ONES_NUM]; 
    bool isRowFound = false;
    U32 nrow = 0;       
#if 0
    U32 r = 1;
    U32 endPoint = L - i - u;
    while (r <= endPoint) //from n row
    {
      for (U32 s = i; s < M; ++s)
      {
        U32 row_ones = 0;
        U32 row_ones_pos[ROW_ONES_NUM];
        U8* row = A[s];
        U32 posIndex = 0;
        for (U32 j = i; j < L - u; ++j) //from take out first i column from I matrix and last u column from U matrix
        {
          if (row[j])
          {
            row_ones++;
            row_ones_pos[posIndex] = j;
            posIndex++;
          }
        }

        if (row_ones == r)
        {
          isRowFound = true;
          nrow  = s;
          memcpy(ones_pos, row_ones_pos, ROW_ONES_NUM);
          break;
        }
      }

      if (isRowFound == true)
      {
        break;
      }
      else
      {
        //Don't have a row with r ones, increase the 1 number of a row to search
        r += 1;
      }
    }
    
    //TODO: Now we must work with this row
    //TODO: if nones == 2 we must choose between the maximum size component (See RFC 5053:5.5.2.2)

    if (isRowFound == false)
    {
      printf("First Phase Decoding Fail!");
      return recovedSym;
    }
#else

    U32 minDegree = m_Degree[i];
    U32 posMinDeg = i;
    U32 r = minDegree;
    if (minDegree == 1)
    {
      r = 1;
    }
    else 
    {     
      if (minDegree == 0)
      {

        if (zerosRowNum > (M - L))
        {
          //free resources
          delete[] c;
          c = NULL;
          
          delete[] d;
          d = NULL;
          
          printf("First Phase Decoding Fail! Zero row Num: %d, M=%d, L=%d\n", zerosRowNum, M, L);
          //return recovedSym;
        }

        U32 lastNonZeroRow = M - zerosRowNum - 1;
        for (U32 zeroIndex = lastNonZeroRow; zeroIndex > i; --zeroIndex)
        {
          if (m_Degree[zeroIndex] != 0)
          {
            lastNonZeroRow = zeroIndex;
            break;
          }
          else
          {
            zerosRowNum++;
          }
        }

        if (i != lastNonZeroRow)
        {
          U8* tmp_row = A[i];
          A[i] = A[lastNonZeroRow];
          A[lastNonZeroRow] = tmp_row;

          zerosRowNum++;

          U32 tmp_encoded = d[i];
          d[i] = d[lastNonZeroRow];
          d[lastNonZeroRow] = tmp_encoded;  

          U32 tmp_degree = m_Degree[i];
          m_Degree[i] = m_Degree[lastNonZeroRow];
          m_Degree[lastNonZeroRow] = tmp_degree;
        }

        minDegree = m_Degree[i];
        posMinDeg = i;
      }

      if (minDegree == 1)
      {
        r = 1;
      }
      else 
      {
        U32 lastNonZeroRow = M - zerosRowNum - 1;
        for (U32 s = i + 1; s <= lastNonZeroRow; ++s)
        {
          if (m_Degree[s] == 1)
          {         
            posMinDeg = s;
            r = 1;
            break;

          }
          else if ((m_Degree[s] > 0) && (m_Degree[s] < minDegree))
          {
            minDegree = m_Degree[s];
            posMinDeg = s;
            r = minDegree;
          }
        }
      }   
    }

    U32 onePosIndex = 0;
    U8 *row = A[posMinDeg];
    nrow = posMinDeg;
    for (U32 j = i; j < L - u; ++j) //from take out first i column from I matrix and last u column from U matrix
    {
      if (row[j])
      {
        ones_pos[onePosIndex] = j;
        onePosIndex++;
      }
    }
#endif

    //Move row to first V row
    if (i != nrow)
    {
      U8* tmp_row = A[i];
      A[i] = A[nrow];
      A[nrow] = tmp_row;

      U32 tmp_encoded = d[i];
      d[i] = d[nrow];
      d[nrow] = tmp_encoded;

      U32 tmp_degree = m_Degree[i];
      m_Degree[i] = m_Degree[nrow];
      m_Degree[nrow] = tmp_degree;
    }

    if (A[i][i] == 0) 
    {
      U32 firstOnePos = ones_pos[0];
      // Change column i and column firstOnePos
      for (U32 rowIndex = 0; rowIndex < M; ++rowIndex)
      {
        U8 tmp_value = A[rowIndex][i];
        A[rowIndex][i] = A[rowIndex][firstOnePos];
        A[rowIndex][firstOnePos] = tmp_value;
      }

      U32 tmp_value = c[i];
      c[i] = c[firstOnePos];
      c[firstOnePos] = tmp_value;
    }

    for (U32 h = 1; h < r; ++h)
    {
      U32 onePos = ones_pos[h];
      U32 changedPos = L - u - 1;
      for (U32 rowIndex = 0; rowIndex < M; ++rowIndex)
      {
        U8 tmp_value = A[rowIndex][changedPos];
        A[rowIndex][changedPos] = A[rowIndex][onePos];
        A[rowIndex][onePos] = tmp_value;
      }

      U32 tmp_value = c[onePos];
      c[onePos] = c[changedPos];
      c[changedPos] = tmp_value;

      u += 1;
    }

    m_Degree[i] -= (r - 1);
    for (U32 k = (i + 1); k < M; ++k)
    {       
      if (A[k][i] == 1)
      {
        m_Degree[k] -= 1; 
      }
      
      for (U32 h = 1; h < r; ++h)
      {
        m_Degree[k] -= A[k][L - u - 1 + h];
      }

      if (A[k][i] == 1)
      {
        for (U32 l = 0; l < L; ++l)
        {
          A[k][l] ^= A[i][l];       
        }

        D[d[k]]->XorData(D[d[i]]);
      }
    }


    printf("\nA matrix is:\n");
    for (U32 step = 0; step < M; ++step)
    {
      for (U32 j = 0; j < L; ++j)
      {
        printf("%d ", A[step][j]);
      }
      printf("--------%d ", m_Degree[step]);
      printf("\n");
    }
    for (U32 i = 0; i < M; ++i)
    {
      U32 length = (*D[i]).GetLen();
      for (U32 j = 0; j < length; ++j)
      {
	printf("%6d", (*D[i]).GetData()[j]);
      }
      if (length==0)
	printf(" 0 ");
    }
    
    i += 1; 
  }
  
  
  /*-----------------------------------Second Phase------------------------------------*/
  U32 rowIndex;
  for (U32 colIndex = i; colIndex < L; ++colIndex)
  {
    rowIndex = colIndex;
    while ((rowIndex < M) && (A[rowIndex][colIndex] == 0))
    {
      rowIndex += 1;
    }
#if 0
    if ((colIndex < (L - i)) && rowIndex == M) 
    {
      printf("Second Phase Decoding Fail!");
      return;
    }
    else if (rowIndex == M)
    {
      continue;
    }
#else
    if (rowIndex == M)
    {
      printf("Second Phase Decoding Fail!");
      //return recovedSym;
    }
#endif

    if (colIndex != rowIndex)
    {
      U8* tmp_row = A[colIndex];
      A[colIndex] = A[rowIndex];
      A[rowIndex] = tmp_row;

      U32 tmp_encoded = d[colIndex];
      d[colIndex] = d[rowIndex];
      d[rowIndex] = tmp_encoded;
    }

    for (U32 k = i; k < M; k++)
    {
      if (k != colIndex)
      {
        if (A[k][colIndex] == 1)
        {
          for (U32 l = i; l < L; l++)
          {
            A[k][l] ^= A[colIndex][l];  
          }

          D[d[k]]->XorData(D[d[colIndex]]);
        }
      }
    }
  }


  printf("\nA matrix before deleting the last M-L rows:\n");
  for (U32 i = 0; i < M; ++i)
  {
    for (U32 j = 0; j < L; ++j)
    {
      printf("%d ", A[i][j]);
    }
    printf("\n");
  }
  
  printf("\nDecoded symbols at first phase:\n");
  for (U32 i = 0; i < M; ++i)
  {
    U32 length = (*D[d[i]]).GetLen();
    for (U32 j = 0; j < length; ++j)
    {
      printf("%6d", (*D[d[i]]).GetData()[j]);
    }
    if (length==0)
      printf(" 0 ");
  }


  for (U32 index = L; index < M; ++index)
  {
    delete[] A[index];
  }

  /*-----------------------------------Third Phase------------------------------------*/
  for (U32 k = 0; k < i; ++k)
  {
    for (U32 s = i; s < L; ++s)
    {
      if (A[k][s] == 1)
      {
        A[k][s] ^= A[s][s];

        D[d[k]]->XorData(D[d[s]]);
      }
    }
  }


  printf("\n");
  printf("\nA Matrix after Gaussian Elimination：\n");
  for (U32 rowIndex = 0; rowIndex < L; rowIndex++)
  {
    for (U32 columnIndex = 0; columnIndex < L; columnIndex++)
    {
      printf("%d ", A[rowIndex][columnIndex]);
    }
    printf("\n");
  }

 
  /*-----------------------------------Fourth Phase------------------------------------*/
  /*
  It is clear at the end of successful decoding that the L symbols D[d[0]], D[d[1]],...,
  D[d[L-1]] are the values of the L symbols C[c[0]], C[c[1]],...,C[c[L-1]]
  */
  for (U32 index = 0; index < L; ++index)
  {
    C[c[index]]->SetData(D[d[index]]);

  }


  printf("\nIntermediate Symbols after decoding LA:\n");
  for (U32 i = 0; i < L; ++i)
  {
    U32 length = (*C[i]).GetLen();
    for (U32 j = 0; j < length; ++j)
    {

      printf("%6d", (*C[i]).GetData()[j]);

    }
  }

}
U16* GetAMatrix(U8** A, U32 M, U32 L)
{
  U16* m_Degree;
  m_Degree = new U16[M];
  memset(m_Degree, 0, M * sizeof(U16));
  printf("\nA matrix row number is：%d，column number：%d\n", M, L);
  for (U32 i = 0; i < M; ++i)
  {
    U16 sum=0;
    for (U32 j = 0; j < L; ++j)
    {
      printf("%d ", A[i][j]);
      if ( A[i][j]==1)
         sum++;
    }
    m_Degree[i]=sum;
    printf("-----%d", m_Degree[i]);
    printf("\n");

    }
  return m_Degree;
}
int main(int argc, char **argv) {
  U32 M=15;
  U32 L=14;  
  U8** A;
  U16* m_Degree; //record the number of 1s in each row
  A = new U8*[M]; //Max number of rows that M can achieve
  //for (U32 i = 0; i < m_L; ++i)
  for (U32 i = 0; i < M; ++i)
  {
    A[i] = new U8[L];
    memset(A[i], 0, L);
  }
  //
    U32 x, y;
  ifstream in("matrix.txt");
  if (!in) {
    cout << "Cannot open file.\n";
    return 0;
  }

  for (x = 0; x < M; x++) {
    for (y = 0; y < L; y++) {
      in >> A[x][y];
    }
  }
  in.close();  
  for (U32 i = 0; i < M; ++i)
  {
    for (U32 j = 0; j < L; ++j)
    {
      A[i][j]=A[i][j]-48;
    }
  }
  
  //
  
  m_Degree= GetAMatrix( A,  M,  L);
  CData** D = new CData*[M]; //Encoded symbols
  CData** C = new CData*[L]; //Encoded symbols
  for (U32 i = 0; i < M; ++i)
  {
    D[i] = new CData();
    if(i<L)
       C[i] = new CData();
    
  } 
  vector<U8> original_data;
  original_data.push_back(0);
  original_data.push_back(0);
  original_data.push_back(0);
  original_data.push_back(0);
  original_data.push_back(0);
  original_data.push_back(0);
  original_data.push_back(0);
  original_data.push_back(0);
  original_data.push_back(0);
  original_data.push_back(0);
  original_data.push_back(81);
  original_data.push_back(255);
  original_data.push_back(74);
  original_data.push_back(236);
  original_data.push_back(220);
  
  //original_data.push_back(45);
  queue<CData*> res;
  U8 *buf = new U8[1];
  U8 dataLen=1;
  for (U32 i = 0; i < M; ++i)
  {
    vector<U8> recover_data ;
    recover_data.push_back(original_data[i]);
    //original_data[i];
    for (U8 i = 0; i < dataLen; ++i)
    {
      buf[i] = recover_data[i];
    }
    res.push(new CData(buf, dataLen));
  }

  U32 EncodedIndexLA= 0;
  while (!res.empty())
  {
  D[EncodedIndexLA]->SetData(res.front());
  res.pop();
  EncodedIndexLA++;
  }
  decode(M, L, A, D,C,m_Degree);  
  return 0;
}