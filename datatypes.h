/* CAUTION: if your data exceeds SIZET_MAX, there would be problem allocating memory*/

#ifndef DATATYPES_INCLUDED

/*datatype for input particle data*/
#ifdef INPUT_REAL8
typedef double IDatReal;
#else
typedef float IDatReal;
#endif

/*datatype for input particle IDs*/
#ifdef INPUT_INT8
typedef long long IDatInt;
#else
#ifdef INPUT_UINT4
typedef unsigned IDatInt;
#else
typedef int IDatInt;
#endif
#endif

/*datatype for internal calculation and output*/
#ifdef HBT_REAL8
typedef double HBTReal;
#else
typedef float HBTReal;
#endif

// the user should ganrantee that HBTInt can at least hold NP_DM (when HBTPID_RANKSTYLE is defined)
#ifdef HBT_INT8
typedef long long HBTInt;  
#define HBTIFMT "%lld"
#else 
typedef int HBTInt;
#define HBTIFMT "%d"
#endif

//auxiliary datatype
typedef HBTReal HBTxyz[3];  //3-d pos/vel data

#if (defined INPUT_REAL8 && defined HBT_REAL8)||(!defined INPUT_REAL8 && !defined HBT_REAL8)
 #define SAME_REALTYPE
#endif

#if (defined INPUT_INT8 && defined HBT_INT8)||(!defined INPUT_INT8 && !defined HBT_INT8)
 #define SAME_INTTYPE
#endif

//constants
#define FRSH_GRPCAT -1
#define FRSH_SUBCAT -2
#define FRSH_SRCCAT -3
//#define FRSH_MBDCAT -4

#define FLAG_LOAD_ID  0b001
#define FLAG_LOAD_POS 0b010
#define FLAG_LOAD_VEL 0b100

#ifndef IniSnap //initial snapshot to start processing
#define IniSnap 0
#endif

#define DATATYPES_INCLUDED
#endif
