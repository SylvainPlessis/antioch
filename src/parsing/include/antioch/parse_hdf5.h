//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-
#ifndef _HDF5_INTERFACE_
#define _HDF5_INTERFACE_

//Antioch
#include "root.h"

//C++
#include <hdf5.h>
#include <string>
#include <vector>
#include <cstdlib>

namespace Antioch
{

namespace root{

class HDF5Loader
{
  public:
     HDF5Loader(){}
     ~HDF5Loader(){}

     hid_t load(const std::string &filename) {file_id = H5Fopen(filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
                                                if(file_id < 0)return file_id;
                                                return H5Gopen(file_id,"/",H5P_DEFAULT);}

    void unload(){if(file_id > 0)H5Fclose(file_id);}

  protected:
    hid_t file_id;                                                

};

class HDF5Reader:public HDF5Loader
{
    public:

      HDF5Reader():myName(""),me(-1),mother(-1){}
      HDF5Reader(const hid_t &above,const hid_t &here, const std::string &name):myName(name),me(here),mother(above){}
      HDF5Reader(const HDF5Reader &rhs);
      ~HDF5Reader(){}

      bool load(const std::string &filename);
      void close(){if(me > 0)H5Gclose(me);me = -1;}

//this attributes
      status queryDoubleAttribute(const std::string &nameAtt, double &att)      const;
      status queryIntAttribute   (const std::string &nameAtt, int &att)         const;
      status queryBoolAttribute  (const std::string &nameAtt, bool &att)        const;
      status queryStringAttribute(const std::string &nameAtt, std::string &att) const;

      status queryVectorDoubleAttribute(const std::string &nameAtt, VectorDouble &att) const;
      status queryVectorIntAttribute   (const std::string &nameAtt, VectorInt &att)    const;
      status queryVectorBoolAttribute  (const std::string &nameAtt, VectorBool &att)   const;
      status queryVectorStringAttribute(const std::string &nameAtt, VectorString &att) const;

      status queryMatrixDoubleAttribute(const std::string &nameAtt, MatrixDouble &att) const;
      status queryMatrixIntAttribute   (const std::string &nameAtt, MatrixInt &att)    const;
      status queryMatrixBoolAttribute  (const std::string &nameAtt, MatrixBool &att)   const;
      status queryMatrixStringAttribute(const std::string &nameAtt, MatrixString &att) const;

//this data
      status queryDoubleData(const std::string &nameAtt, double &att)      const;
      status queryIntData   (const std::string &nameAtt, int &att)         const;
      status queryBoolData  (const std::string &nameAtt, bool &att)        const;
      status queryStringData(const std::string &nameAtt, std::string &att) const;

      status queryVectorDoubleData(const std::string &nameAtt, VectorDouble &att) const;
      status queryVectorIntData   (const std::string &nameAtt, VectorInt &att)    const;
      status queryVectorBoolData  (const std::string &nameAtt, VectorBool &att)   const;
      status queryVectorStringData(const std::string &nameAtt, VectorString &att) const;

      status queryMatrixDoubleData(const std::string &nameAtt, MatrixDouble &att) const;
      status queryMatrixIntData   (const std::string &nameAtt, MatrixInt &att)    const;
      status queryMatrixBoolData  (const std::string &nameAtt, MatrixBool &att)   const;
      status queryMatrixStringData(const std::string &nameAtt, MatrixString &att) const;

//this data's attribute
      status queryDoubleDataAttribute(const std::string &nameData, const std::string &nameAtt, double &att)      const;
      status queryIntDataAttribute   (const std::string &nameData, const std::string &nameAtt, int &att)         const;
      status queryBoolDataAttribute  (const std::string &nameData, const std::string &nameAtt, bool &att)        const;
      status queryStringDataAttribute(const std::string &nameData, const std::string &nameAtt, std::string &att) const;

      status queryVectorDoubleDataAttribute(const std::string &nameData, const std::string &nameAtt, VectorDouble &att) const;
      status queryVectorIntDataAttribute   (const std::string &nameData, const std::string &nameAtt, VectorInt &att)    const;
      status queryVectorBoolDataAttribute  (const std::string &nameData, const std::string &nameAtt, VectorBool &att)   const;
      status queryVectorStringDataAttribute(const std::string &nameData, const std::string &nameAtt, VectorString &att) const;

      status queryMatrixDoubleDataAttribute(const std::string &nameData, const std::string &nameAtt, MatrixDouble &att) const;
      status queryMatrixIntDataAttribute   (const std::string &nameData, const std::string &nameAtt, MatrixInt &att)    const;
      status queryMatrixBoolDataAttribute  (const std::string &nameData, const std::string &nameAtt, MatrixBool &att)   const;
      status queryMatrixStringDataAttribute(const std::string &nameData, const std::string &nameAtt, MatrixString &att) const;

      HDF5Reader& operator= (const HDF5Reader &rhs);

      const hid_t getMe()       const {return me;}
      const hid_t getMom()      const {return mother;}
      const std::string getId() const {return myName;}

      const HDF5Reader firstChild()  const;
      const HDF5Reader nextSibling() const;
    private:
//data
      template<typename readType>
      status queryData(const std::string &nameAtt, readType *&value, unsigned int &sizex, unsigned int &sizey, bool &ischar) const;

      template<typename outType, typename readType>
      status queryScalarData(const std::string &nameAtt, outType &value) const;
      template<typename outType, typename readType>
      status queryVectorData(const std::string &nameAtt, outType &value) const;
      template<typename outType, typename readType>
      status queryMatrixData(const std::string &nameAtt, outType &value) const;

//attributes
      template<typename readType>
      status queryAttribute(const std::string &nameAtt, readType *&value, unsigned int &sizex, unsigned int &sizey, bool &ischar, const hid_t &where) const;

      template<typename outType, typename readType>
      status queryScalarAttribute(const std::string &nameAtt, outType &value, hid_t where = 0) const;
      template<typename outType, typename readType>
      status queryVectorAttribute(const std::string &nameAtt, outType &value, hid_t where = 0) const;
      template<typename outType, typename readType>
      status queryMatrixAttribute(const std::string &nameAtt, outType &value, hid_t where = 0) const;

//data attributes
      template<typename TypeIn, typename TypeOut>
      status queryScalarDataAttribute(const std::string &nameData, const std::string &nameAtt, TypeIn &value) const;
      template<typename TypeIn, typename TypeOut>
      status queryVectorDataAttribute(const std::string &nameData, const std::string &nameAtt, TypeIn &value) const;
      template<typename TypeIn, typename TypeOut>
      status queryMatrixDataAttribute(const std::string &nameData, const std::string &nameAtt, TypeIn &value) const;


      void deleteifchar(void *reader) const;

      std::string myName;
      hid_t me;
      hid_t mother;
};

class HDF5Writer: public HDF5Loader
{
    public:

      HDF5Writer(const hid_t &above,const hid_t &here,const std::string &name):me(here),mother(above),myName(name){}
      HDF5Writer(const HDF5Writer &rhs);
      ~HDF5Writer(){if(me > 0)H5Gclose(me);}

      bool load(const std::string &filename);

//this attributes
      status writeDoubleAttribute(const std::string &nameAtt, const double &att)      {return writeScalarAttribute<double,double>(nameAtt,att,DOUBLE,me);}
      status writeIntAttribute   (const std::string &nameAtt, const int &att)         {return writeScalarAttribute<int,int>   (nameAtt,att,INT,me);}
      status writeBoolAttribute  (const std::string &nameAtt, const bool &att)        {return writeScalarAttribute<bool,int>   (nameAtt,att,INT,me);}
      status writeStringAttribute(const std::string &nameAtt, const std::string &att) {return writeScalarAttribute<std::string,std::string>(nameAtt,att,STRING,me);}

      status writeVectorDoubleAttribute(const std::string &nameAtt, const VectorDouble &att) {return writeVectorAttribute<VectorDouble,double>(nameAtt,att,DOUBLE,me);}
      status writeVectorIntAttribute   (const std::string &nameAtt, const VectorInt &att)    {return writeVectorAttribute<VectorInt,int>(nameAtt,att,INT,me);}
      status writeVectorBoolAttribute  (const std::string &nameAtt, const VectorBool &att)   {return writeVectorAttribute<VectorBool,int>(nameAtt,att,INT,me);}
      status writeVectorStringAttribute(const std::string &nameAtt, const VectorString &att) {return writeVectorAttribute<VectorString,std::string>(nameAtt,att,STRING,me);}

      status writeMatrixDoubleAttribute(const std::string &nameAtt, const MatrixDouble &att){return writeMatrixAttribute<MatrixDouble,double>(nameAtt,att,DOUBLE,me);}
      status writeMatrixIntAttribute   (const std::string &nameAtt, const MatrixInt &att)   {return writeMatrixAttribute<MatrixInt,int>(nameAtt,att,INT,me);}
      status writeMatrixBoolAttribute  (const std::string &nameAtt, const MatrixBool &att)  {return writeMatrixAttribute<MatrixBool,int>(nameAtt,att,INT,me);}
      status writeMatrixStringAttribute(const std::string &nameAtt, const MatrixString &att){return writeMatrixAttribute<MatrixString,std::string>(nameAtt,att,STRING,me);}

//this data
      status writeDoubleData(const std::string &nameAtt, const double &att)      {return writeScalarData<double,double>(nameAtt,att,DOUBLE,me);}
      status writeIntData   (const std::string &nameAtt, const int &att)         {return writeScalarData<int,int>(nameAtt,att,INT,me);}
      status writeBoolData  (const std::string &nameAtt, const bool &att)        {return writeScalarData<bool,int>(nameAtt,att,INT,me);}
      status writeStringData(const std::string &nameAtt, const std::string &att) {return writeScalarData<std::string,std::string>(nameAtt,att,STRING,me);}

      status writeVectorDoubleData(const std::string &nameAtt, const VectorDouble &att) {return writeVectorData<VectorDouble,double>(nameAtt,att,DOUBLE,me);}
      status writeVectorIntData   (const std::string &nameAtt, const VectorInt &att)    {return writeVectorData<VectorInt,int>(nameAtt,att,INT,me);}
      status writeVectorBoolData  (const std::string &nameAtt, const VectorBool &att)   {return writeVectorData<VectorBool,int>(nameAtt,att,INT,me);}
      status writeVectorStringData(const std::string &nameAtt, const VectorString &att) {return writeVectorData<VectorString,std::string>(nameAtt,att,STRING,me);}

      status writeMatrixDoubleData(const std::string &nameAtt, const MatrixDouble &att) {return writeMatrixData<MatrixDouble,double>(nameAtt,att,DOUBLE,me);}
      status writeMatrixIntData   (const std::string &nameAtt, const MatrixInt &att)    {return writeMatrixData<MatrixInt,int>(nameAtt,att,INT,me);}
      status writeMatrixBoolData  (const std::string &nameAtt, const MatrixBool &att)   {return writeMatrixData<MatrixBool,int>(nameAtt,att,INT,me);}
      status writeMatrixStringData(const std::string &nameAtt, const MatrixString &att) {return writeMatrixData<MatrixString,std::string>(nameAtt,att,STRING,me);}

//this data's attribute
      status writeDoubleDataAttribute(const std::string &nameData, const std::string &nameAtt, const double &att);
      status writeIntDataAttribute   (const std::string &nameData, const std::string &nameAtt, const int &att);
      status writeBoolDataAttribute  (const std::string &nameData, const std::string &nameAtt, const bool &att);
      status writeStringDataAttribute(const std::string &nameData, const std::string &nameAtt, const std::string &att);

      status writeVectorDoubleDataAttribute(const std::string &nameData, const std::string &nameAtt, const VectorDouble &att);
      status writeVectorIntDataAttribute   (const std::string &nameData, const std::string &nameAtt, const VectorInt &att);
      status writeVectorBoolDataAttribute  (const std::string &nameData, const std::string &nameAtt, const VectorBool &att);
      status writeVectorStringDataAttribute(const std::string &nameData, const std::string &nameAtt, const VectorString &att);

      status writeMatrixDoubleDataAttribute(const std::string &nameData, const std::string &nameAtt, const MatrixDouble &att);
      status writeMatrixIntDataAttribute   (const std::string &nameData, const std::string &nameAtt, const MatrixInt &att);
      status writeMatrixBoolDataAttribute  (const std::string &nameData, const std::string &nameAtt, const MatrixBool &att);
      status writeMatrixStringDataAttribute(const std::string &nameData, const std::string &nameAtt, const MatrixString &att);

      HDF5Writer& operator= (const HDF5Writer &rhs);

      const hid_t getMe()       const {return me;}
      const hid_t getMom()      const {return mother;}
      const std::string getId() const {return myName;}

      HDF5Writer addChild(const std::string &child);

    private:
      hid_t me;
      hid_t mother;
      std::string myName;

      enum typeenum{DOUBLE = 0,INT,STRING};

//this attribute
      template<typename TypeIn, typename TypeOut>
      status writeScalarAttribute(const std::string &nameAtt, const TypeIn &att, const typeenum &type, const hid_t &mother);
      template<typename TypeIn, typename TypeOut>
      status writeVectorAttribute(const std::string &nameAtt, const TypeIn &att, const typeenum &type, const hid_t &mother);
      template<typename TypeIn, typename TypeOut>
      status writeMatrixAttribute(const std::string &nameAtt, const TypeIn &att, const typeenum &type, const hid_t &mother);

      template<typename TypeAtt>
      status writeAttribute(const std::string &nameAtt, const TypeAtt *att, unsigned int sizex, unsigned int sizey, const hid_t &where, const typeenum &type);

//data
      template<typename TypeDat>
      status writeData(const std::string &nameData, const TypeDat *data, unsigned int sizex, unsigned int sizey, const hid_t &where, const typeenum &type);

      template<typename TypeIn, typename TypeOut>
      status writeScalarData(const std::string &nameData, const TypeIn &dat, const typeenum &type, const hid_t &mother);
      template<typename TypeIn, typename TypeOut>
      status writeVectorData(const std::string &nameData, const TypeIn &dat, const typeenum &type, const hid_t &mother);
      template<typename TypeIn, typename TypeOut>
      status writeMatrixData(const std::string &nameData, const TypeIn &dat, const typeenum &type, const hid_t &mother);

//data attribute
      template<typename TypeIn, typename TypeOut>
      status writeScalarDataAttribute(const std::string &nameData, const std::string &nameAtt, const TypeIn &att, const typeenum &type);
      template<typename TypeIn, typename TypeOut>
      status writeVectorDataAttribute(const std::string &nameData, const std::string &nameAtt, const TypeIn &att, const typeenum &type);
      template<typename TypeIn, typename TypeOut>
      status writeMatrixDataAttribute(const std::string &nameData, const std::string &nameAtt, const TypeIn &att, const typeenum &type);
};

inline
HDF5Writer::HDF5Writer(const HDF5Writer &rhs)
{
  *this = rhs;
}

inline
HDF5Writer &HDF5Writer::operator=(const HDF5Writer &rhs)
{
  if(this == &rhs)return *this;
  me = rhs.getMe();
  mother = rhs.getMom();
  myName = rhs.getId();
  return *this;
}

inline
bool HDF5Writer::load(const std::string &filename)
{
  me = HDF5Loader::load(filename);
  mother = file_id;
  myName = "/";
  return (me > 0);
}

inline
HDF5Writer HDF5Writer::addChild(const std::string &child)
{
  hid_t c = H5Gopen2(me,child.c_str(),H5P_DEFAULT);

  return HDF5Writer(me,c,child);
}

//attribute
template<typename TypeIn, typename TypeOut>
inline
status HDF5Writer::writeMatrixAttribute(const std::string &nameAtt, const TypeIn &att, const typeenum &type, const hid_t &where)
{
   unsigned int sizex(att.size()),sizey(att.front().size());
   TypeOut toWrite[sizex * sizey];
   for(unsigned int i = 0; i < sizex; i++)
   {
     for(unsigned int j = 0; j < sizey; j++)
     {
       toWrite[j + i*sizex] = att[i][j];
     }
   }

   return writeAttribute<TypeOut>(nameAtt,toWrite,sizex,sizey,where,type);
}

template<typename TypeIn, typename TypeOut >
inline
status HDF5Writer::writeVectorAttribute(const std::string &nameAtt, const TypeIn &att, const typeenum &type, const hid_t &where)
{
   unsigned int sizex(att.size()),sizey(0);
   TypeOut toWrite[sizex];
   for(unsigned int i = 0; i < sizex; i++)
   {
     toWrite[i] = att[i];
   }

   return writeAttribute<TypeOut>(nameAtt,toWrite,sizex,sizey,where,type);

}

template<typename TypeIn, typename TypeOut>
inline
status HDF5Writer::writeScalarAttribute(const std::string &nameAtt, const TypeIn &att, const typeenum &type, const hid_t &where)
{
   unsigned int sizex(1),sizey(0);
   TypeOut toWrite[sizex];
   toWrite[0] = att;

   return writeAttribute<TypeOut>(nameAtt,toWrite,sizex,sizey,where,type);
}

template <typename TypeAtt>
inline
status HDF5Writer::writeAttribute(const std::string &nameAtt, const TypeAtt *att, unsigned int sizex, unsigned int sizey, const hid_t &where,
                                  const typeenum &type)
{
  hid_t type_mem(-1);
  herr_t stat;
  switch(type)
  {
   case DOUBLE:
     type_mem = H5T_NATIVE_DOUBLE;
     break;
   case INT:
     type_mem = H5T_NATIVE_INT;
     break;
   case STRING:
     type_mem = H5Tcopy(H5T_C_S1);
     stat = H5Tset_size(type_mem,H5T_VARIABLE);
     break;
    default:
     return ERROR;
     break;
  }
  hsize_t Dim(1);
  if(sizey != 0)Dim = 2;
  hsize_t Dims[Dim];
  Dims[0] = sizex;
  if(sizey != 0)Dims[1] = sizey;

  hid_t space_id = H5Screate_simple(Dim,Dims,NULL);
  hid_t att_id = H5Acreate2(where,nameAtt.c_str(),type_mem,space_id,H5P_DEFAULT,H5P_DEFAULT);
  stat = H5Awrite(att_id,type_mem,att);
  stat = H5Aclose(att_id);
  stat = H5Sclose(space_id);
  if(type == STRING)H5Tclose(type_mem);

  return SUCCESS;
}

//data
template<typename TypeIn, typename TypeOut>
inline
status HDF5Writer::writeMatrixData(const std::string &nameData, const TypeIn &data, const typeenum &type, const hid_t &where)
{
   unsigned int sizex(data.size()),sizey(data.front().size());
   TypeOut toWrite[sizex * sizey];
   for(unsigned int i = 0; i < sizex; i++)
   {
     for(unsigned int j = 0; j < sizey; j++)
     {
       toWrite[i*sizex + j] = data[i][j];
     }
   }

   return writeData<TypeOut>(nameData,toWrite,sizex,sizey,where,type);
}

template<typename TypeIn, typename TypeOut>
inline
status HDF5Writer::writeVectorData(const std::string &nameData, const TypeIn &data, const typeenum &type, const hid_t &where)
{
   unsigned int sizex(data.size()),sizey(0);
   TypeOut toWrite[sizex];
   for(unsigned int i = 0; i < sizex; i++)
   {
     toWrite[i] = data[i];
   }

   return writeData<TypeOut>(nameData,toWrite,sizex,sizey,where,type);

}

template<typename TypeIn, typename TypeOut>
inline
status HDF5Writer::writeScalarData(const std::string &nameData, const TypeIn &data, const typeenum &type, const hid_t &where)
{
   unsigned int sizex(1),sizey(0);
   TypeOut toWrite[sizex];
   toWrite[0] = data;

   return writeData<TypeOut>(nameData,toWrite,sizex,sizey,where,type);

}

template <typename TypeData>
inline
status HDF5Writer::writeData(const std::string &nameData, const TypeData *data, unsigned int sizex, unsigned int sizey, const hid_t &where,
                                  const typeenum &type)
{
  hid_t type_mem(-1);
  herr_t stat;
  switch(type)
  {
   case DOUBLE:
     type_mem = H5T_NATIVE_DOUBLE;
     break;
   case INT:
     type_mem = H5T_NATIVE_INT;
     break;
   case STRING:
     type_mem = H5Tcopy(H5T_C_S1);
     stat = H5Tset_size(type_mem,H5T_VARIABLE);
     break;
    default:
     return ERROR;
     break;
  }
  hsize_t Dim(1);
  if(sizey != 0)Dim = 2;
  hsize_t Dims[Dim];
  Dims[0] = sizex;
  if(sizey != 0)Dims[1] = sizey;

  hid_t space_id = H5Screate_simple(Dim,Dims,NULL);
  hid_t data_id = H5Dcreate2(where,nameData.c_str(),type_mem,space_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  stat = H5Dwrite(data_id,type_mem,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  stat = H5Dclose(data_id);
  stat = H5Sclose(space_id);

  return SUCCESS;
}

//data's attribute
template <typename TypeIn, typename TypeOut>
inline
status HDF5Writer::writeScalarDataAttribute(const std::string &nameData, const std::string &nameAtt, const TypeIn &att, const typeenum &type)
{
  hid_t data_id = H5Dopen2(me,nameData.c_str(),H5P_DEFAULT);
  if(data_id < 0)return NOT_FOUND;
  status stat = writeScalarAttribute<TypeIn,TypeOut>(nameAtt,att,type,data_id);
  H5Dclose(data_id);
  return stat;
}

template <typename TypeIn, typename TypeOut>
inline
status HDF5Writer::writeVectorDataAttribute(const std::string &nameData, const std::string &nameAtt, const TypeIn &att, const typeenum &type)
{
  hid_t data_id = H5Dopen2(me,nameData.c_str(),H5P_DEFAULT);
  if(data_id < 0)return NOT_FOUND;
  status stat = writeVectorAttribute<TypeIn,TypeOut>(nameAtt,att,type,data_id);
  H5Dclose(data_id);
  return stat;
}

template <typename TypeIn, typename TypeOut>
inline
status HDF5Writer::writeMatrixDataAttribute(const std::string &nameData, const std::string &nameAtt, const TypeIn &att, const typeenum &type)
{
  hid_t data_id = H5Dopen2(me,nameData.c_str(),H5P_DEFAULT);
  if(data_id < 0)return NOT_FOUND;
  status stat = writeMatrixAttribute<TypeIn,TypeOut>(nameAtt,att,type,data_id);
  H5Dclose(data_id);
  return stat;
}

inline
status HDF5Writer::writeDoubleDataAttribute(const std::string &nameData, const std::string &nameAtt, const double &att)
{
  return writeScalarDataAttribute<double,double>(nameData,nameAtt,att,DOUBLE);
}
inline
status HDF5Writer::writeIntDataAttribute(const std::string &nameData, const std::string &nameAtt, const int &att)
{
  return writeScalarDataAttribute<int,int>(nameData,nameAtt,att,INT);
}
inline
status HDF5Writer::writeBoolDataAttribute(const std::string &nameData, const std::string &nameAtt, const bool &att)
{
  return writeScalarDataAttribute<bool,int>(nameData,nameAtt,att,INT);
}
inline
status HDF5Writer::writeStringDataAttribute(const std::string &nameData, const std::string &nameAtt, const std::string &att)
{
  return writeScalarDataAttribute<std::string,std::string>(nameData,nameAtt,att,STRING);
}

inline
status HDF5Writer::writeVectorDoubleDataAttribute(const std::string &nameData, const std::string &nameAtt, const VectorDouble &att)
{
  return writeVectorDataAttribute<VectorDouble,double>(nameData,nameAtt,att,DOUBLE);
}
inline
status HDF5Writer::writeVectorIntDataAttribute(const std::string &nameData, const std::string &nameAtt, const VectorInt &att)
{
  return writeVectorDataAttribute<VectorInt,int>(nameData,nameAtt,att,INT);
}
inline
status HDF5Writer::writeVectorBoolDataAttribute  (const std::string &nameData, const std::string &nameAtt, const VectorBool &att)
{
  return writeVectorDataAttribute<VectorBool,int>(nameData,nameAtt,att,INT);
}
inline
status HDF5Writer::writeVectorStringDataAttribute(const std::string &nameData, const std::string &nameAtt, const VectorString &att)
{
  return writeVectorDataAttribute<VectorString,std::string>(nameData,nameAtt,att,STRING);
}

inline
status HDF5Writer::writeMatrixDoubleDataAttribute(const std::string &nameData, const std::string &nameAtt, const MatrixDouble &att)
{
  return writeMatrixDataAttribute<MatrixDouble,double>(nameData,nameAtt,att,DOUBLE);
}
inline
status HDF5Writer::writeMatrixIntDataAttribute(const std::string &nameData, const std::string &nameAtt, const MatrixInt &att)
{
  return writeMatrixDataAttribute<MatrixInt,int>(nameData,nameAtt,att,INT);
}
inline
status HDF5Writer::writeMatrixBoolDataAttribute(const std::string &nameData, const std::string &nameAtt, const MatrixBool &att)
{
  return writeMatrixDataAttribute<MatrixBool,int>(nameData,nameAtt,att,INT);
}
inline
status HDF5Writer::writeMatrixStringDataAttribute(const std::string &nameData, const std::string &nameAtt, const MatrixString &att)
{
  return writeMatrixDataAttribute<MatrixString,std::string>(nameData,nameAtt,att,STRING);
}

////////////////////////////////
//
// reader's methods
// 
///////////////////////////////
inline
bool HDF5Reader::load(const std::string &filename)
{
  me = HDF5Loader::load(filename);
  mother = file_id;
  myName = "/";
  return (me > 0);
}

inline
const HDF5Reader HDF5Reader::firstChild() const
{
   H5G_info_t info;
   herr_t hstatus = H5Gget_info(me,&info);
   if(hstatus < 0)throw ERROR;
   for(unsigned int i = 0; i < info.nlinks; i++)
   {
      char nameDataC[1024];
      H5O_info_t objectInfo;
      hstatus = H5Oget_info_by_idx(me,".", 
                                  H5_INDEX_NAME,H5_ITER_NATIVE,
                                  i, &objectInfo,H5P_DEFAULT); // to get the type of object
      H5Lget_name_by_idx(me,".", 
                         H5_INDEX_NAME,H5_ITER_NATIVE,
                         i, nameDataC,1024,H5P_DEFAULT); // to get the name
      std::string nameObject(nameDataC);
      if(objectInfo.type == H5O_TYPE_GROUP)
      {
         hid_t group_id = H5Gopen2(me,nameObject.c_str(),H5P_DEFAULT);
         HDF5Reader child(me,group_id,nameObject);
         return child;
      }
   }
   throw NOT_FOUND;
   return HDF5Reader();
}

inline
const HDF5Reader HDF5Reader::nextSibling() const
{
   H5G_info_t info;
   herr_t hstatus = H5Gget_info(mother,&info);
   if(hstatus < 0)throw ERROR;
   bool yet(false);
   for(unsigned int i = 0; i < info.nlinks; i++)
   {
      char nameDataC[1024];
      H5O_info_t objectInfo;
      hstatus = H5Oget_info_by_idx(mother,".", 
                                  H5_INDEX_NAME,H5_ITER_NATIVE ,
                                  i, &objectInfo,H5P_DEFAULT); // to get the type of object
      H5Lget_name_by_idx(mother,".", 
                         H5_INDEX_NAME,H5_ITER_NATIVE ,
                         i, nameDataC,1024,H5P_DEFAULT); // to get the name
      std::string nameObject(nameDataC);
      if(nameObject == myName)
      {
        yet = true;
        continue;
      }
      if(!yet)continue;
      if(objectInfo.type == H5O_TYPE_GROUP)
      {
         hid_t group_id = H5Gopen2(mother,nameObject.c_str(),H5P_DEFAULT);
         HDF5Reader sibling(mother,group_id,nameObject);
         return sibling;
      }
   }
   throw NOT_FOUND;
   return HDF5Reader();
}

inline
status HDF5Reader::queryDoubleAttribute(const std::string &nameAtt, double &att) const
{
  return queryScalarAttribute<double,double>(nameAtt,att);
}

inline
status HDF5Reader::queryIntAttribute(const std::string &nameAtt, int &att) const
{
  return queryScalarAttribute<int,int>(nameAtt,att);
}

inline
status HDF5Reader::queryBoolAttribute(const std::string &nameAtt, bool &att) const
{
  return queryScalarAttribute<bool,int>(nameAtt,att);
}

inline
status HDF5Reader::queryStringAttribute(const std::string &nameAtt, std::string &att) const
{
  return queryScalarAttribute<std::string,char*>(nameAtt,att);
}

inline
status HDF5Reader::queryVectorDoubleAttribute(const std::string &nameAtt, VectorDouble &att) const
{
  return queryVectorAttribute<VectorDouble,double>(nameAtt,att);
}

inline
status HDF5Reader::queryVectorIntAttribute(const std::string &nameAtt, VectorInt &att) const
{
  return queryVectorAttribute<VectorInt,int>(nameAtt,att);
}

inline
status HDF5Reader::queryVectorBoolAttribute(const std::string &nameAtt, VectorBool &att) const
{
  return queryVectorAttribute<VectorBool,int>(nameAtt,att);
}

inline
status HDF5Reader::queryVectorStringAttribute(const std::string &nameAtt, VectorString &att) const
{
  return queryVectorAttribute<VectorString,char*>(nameAtt,att);
}

inline
status HDF5Reader::queryMatrixDoubleAttribute(const std::string &nameAtt, MatrixDouble &att) const
{
  return queryMatrixAttribute<MatrixDouble,double>(nameAtt,att);
}

inline
status HDF5Reader::queryMatrixIntAttribute(const std::string &nameAtt, MatrixInt &att) const
{
  return queryMatrixAttribute<MatrixInt,int>(nameAtt,att);
}

inline
status HDF5Reader::queryMatrixBoolAttribute(const std::string &nameAtt, MatrixBool &att) const
{
  return queryMatrixAttribute<MatrixBool,int>(nameAtt,att);
}

inline
status HDF5Reader::queryMatrixStringAttribute(const std::string &nameAtt, MatrixString &att) const
{
  return queryMatrixAttribute<MatrixString,char*>(nameAtt,att);
}

//this data
inline
status HDF5Reader::queryDoubleData(const std::string &nameAtt, double &att) const
{
  return queryScalarData<double,double>(nameAtt,att);
}

inline
status HDF5Reader::queryIntData(const std::string &nameAtt, int &att) const
{
  return queryScalarData<int,int>(nameAtt,att);
}

inline
status HDF5Reader::queryBoolData(const std::string &nameAtt, bool &att) const
{
  return queryScalarData<bool,int>(nameAtt,att);
}

inline
status HDF5Reader::queryStringData(const std::string &nameAtt, std::string &att) const
{
  return queryScalarData<std::string,char*>(nameAtt,att);
}


inline
status HDF5Reader::queryVectorDoubleData(const std::string &nameAtt, VectorDouble &att) const
{
  return queryVectorData<VectorDouble,double>(nameAtt,att);
}

inline
status HDF5Reader::queryVectorIntData(const std::string &nameAtt, VectorInt &att) const
{
  return queryVectorData<VectorInt,int>(nameAtt,att);
}

inline
status HDF5Reader::queryVectorBoolData(const std::string &nameAtt, VectorBool &att) const
{
  return queryVectorData<VectorBool,int>(nameAtt,att);
}

inline
status HDF5Reader::queryVectorStringData(const std::string &nameAtt, VectorString &att) const
{
  return queryVectorData<VectorString,char*>(nameAtt,att);
}

inline
status HDF5Reader::queryMatrixDoubleData(const std::string &nameAtt, MatrixDouble &att) const
{
  return queryMatrixData<MatrixDouble,double>(nameAtt,att);
}

inline
status HDF5Reader::queryMatrixIntData(const std::string &nameAtt, MatrixInt &att) const
{
  return queryMatrixData<MatrixInt,int>(nameAtt,att);
}

inline
status HDF5Reader::queryMatrixBoolData(const std::string &nameAtt, MatrixBool &att) const
{
  return queryMatrixData<MatrixBool,int>(nameAtt,att);
}

inline
status HDF5Reader::queryMatrixStringData(const std::string &nameAtt, MatrixString &att)  const
{
  return queryMatrixData<MatrixString,char*> (nameAtt,att);
}

//this data's attribute
template <typename TypeIn, typename TypeOut>
inline
status HDF5Reader::queryScalarDataAttribute(const std::string &nameData, const std::string &nameAtt, TypeIn &att) const
{
  hid_t data_id = H5Dopen2(me,nameData.c_str(),H5P_DEFAULT);
  if(data_id < 0)return NOT_FOUND;
  status stat = queryScalarAttribute<TypeIn,TypeOut>(nameAtt,att,data_id);
  H5Dclose(data_id);
  return stat;
}

template <typename TypeIn, typename TypeOut>
inline
status HDF5Reader::queryVectorDataAttribute(const std::string &nameData, const std::string &nameAtt, TypeIn &att) const
{
  hid_t data_id = H5Dopen2(me,nameData.c_str(),H5P_DEFAULT);
  if(data_id < 0)return NOT_FOUND;
  status stat = queryVectorAttribute<TypeIn,TypeOut>(nameAtt,att,data_id);
  H5Dclose(data_id);
  return stat;
}

template <typename TypeIn, typename TypeOut>
inline
status HDF5Reader::queryMatrixDataAttribute(const std::string &nameData, const std::string &nameAtt, TypeIn &att) const
{
  hid_t data_id = H5Dopen2(me,nameData.c_str(),H5P_DEFAULT);
  if(data_id < 0)return NOT_FOUND;
  status stat = queryMatrixAttribute<TypeIn,TypeOut>(nameAtt,att,data_id);
  H5Dclose(data_id);
  return stat;
}

inline
status HDF5Reader::queryDoubleDataAttribute(const std::string &nameData, const std::string &nameAtt, double &att) const
{
  return queryScalarDataAttribute<double,double>(nameData,nameAtt,att);
}

inline
status HDF5Reader::queryIntDataAttribute(const std::string &nameData, const std::string &nameAtt, int &att) const
{
  return queryScalarDataAttribute<int,int>(nameData,nameAtt,att);
}

inline
status HDF5Reader::queryBoolDataAttribute(const std::string &nameData, const std::string &nameAtt, bool &att) const
{
  return queryScalarDataAttribute<bool,int>(nameData,nameAtt,att);
}

inline
status HDF5Reader::queryStringDataAttribute(const std::string &nameData, const std::string &nameAtt, std::string &att) const
{
  return queryScalarDataAttribute<std::string,char*>(nameData,nameAtt,att);
}


inline
status HDF5Reader::queryVectorDoubleDataAttribute(const std::string &nameData, const std::string &nameAtt, VectorDouble &att) const
{
  return queryVectorDataAttribute<VectorDouble,double>(nameData,nameAtt,att);
}

inline
status HDF5Reader::queryVectorIntDataAttribute(const std::string &nameData, const std::string &nameAtt, VectorInt &att) const
{
  return queryVectorDataAttribute<VectorInt,int>(nameData,nameAtt,att);
}

inline
status HDF5Reader::queryVectorBoolDataAttribute(const std::string &nameData, const std::string &nameAtt, VectorBool &att) const
{
  return queryVectorDataAttribute<VectorBool,int>(nameData,nameAtt,att);
}

inline
status HDF5Reader::queryVectorStringDataAttribute(const std::string &nameData, const std::string &nameAtt, VectorString &att) const
{
  return queryVectorDataAttribute<VectorString,char*>(nameData,nameAtt,att);
}


inline
status HDF5Reader::queryMatrixDoubleDataAttribute(const std::string &nameData, const std::string &nameAtt, MatrixDouble &att) const
{
  return queryMatrixDataAttribute<MatrixDouble,double>(nameData,nameAtt,att);
}

inline
status HDF5Reader::queryMatrixIntDataAttribute(const std::string &nameData, const std::string &nameAtt, MatrixInt &att) const
{
  return queryMatrixDataAttribute<MatrixInt,int>(nameData,nameAtt,att);
}

inline
status HDF5Reader::queryMatrixBoolDataAttribute(const std::string &nameData, const std::string &nameAtt, MatrixBool &att) const
{
  return queryMatrixDataAttribute<MatrixBool,int>(nameData,nameAtt,att);
}

inline
status HDF5Reader::queryMatrixStringDataAttribute(const std::string &nameData, const std::string &nameAtt, MatrixString &att) const
{
  return queryMatrixDataAttribute<MatrixString,char*>(nameData,nameAtt,att);
}


inline
void HDF5Reader::deleteifchar(void *reader) const
{
  char *charprt = static_cast<char*>(*static_cast<char**>(reader));
  free(charprt);
}

inline
HDF5Reader &HDF5Reader::operator=(const HDF5Reader &rhs)
{
  if(this == &rhs)return *this;
  me = rhs.getMe();
  mother = rhs.getMom();
  myName = rhs.getId();
  return *this;
}

inline
HDF5Reader::HDF5Reader(const HDF5Reader &rhs)
{
  *this = rhs;
}

template <typename readType>
inline
status HDF5Reader::queryData(const std::string &nameData, readType *&reader, 
                                              unsigned int &sizex, unsigned int &sizey, 
                                              bool &ischar) const
{
   hid_t data_id       = H5Dopen2(me,nameData.c_str(),H5P_DEFAULT);
   if(data_id <= 0)return ERROR;
   hid_t data_space_id = H5Dget_space(data_id); // set space to obtain sizes
   hid_t data_type_id  = H5Dget_type(data_id);

   int nrank = H5Sget_simple_extent_ndims(data_space_id); // number of rows (values per parameters)
   if(nrank > 2)return ERROR;
   hsize_t dims[nrank]; // nrankData, dimension of tensor
   herr_t hstatus = H5Sget_simple_extent_dims(data_space_id,dims,NULL); // fills with each tensor dimension's dimension
   if(hstatus < 0)return ERROR;
   hid_t data_space_mem_id = H5Screate_simple(nrank,dims,NULL);
   unsigned int fullSize(dims[0]),width(0);
   if(nrank == 2)
   {
     fullSize *= dims[1];
     width = dims[1];
   }
   sizex = dims[0];
   sizey = width;
   reader = new readType[fullSize];

//dims[1]*y + x
//in case of string
   hid_t string_type = H5Tcopy(H5T_C_S1);
   if(H5Tget_class(data_type_id) == H5T_STRING)
   {
     ischar = true;
     if(H5Tis_variable_str(data_type_id))
     {
       string_type = H5Tget_native_type(data_type_id, H5T_DIR_ASCEND);
     }else
     {
       hstatus = H5Tset_size(string_type,H5Tget_size(data_type_id));
       if(hstatus < 0)return ERROR;
     }
   }

   hstatus = H5Dread(data_id,data_type_id,data_space_id,data_space_mem_id,H5P_DEFAULT,reader);

   if(hstatus < 0)return ERROR;

   if(H5Tget_class(data_type_id) == H5T_STRING)H5Tclose(string_type);

//close
   H5Tclose(data_type_id);
   H5Dclose(data_id);
   H5Sclose(data_space_id);
   H5Sclose(data_space_mem_id);

   return SUCCESS;

}

template<typename outType, typename readType>
inline
status HDF5Reader::queryScalarData(const std::string &nameAtt, outType &value) const
{
   readType *reader(NULL);
   unsigned int sizex(0),sizey(0);
   bool ischar(false);
   status out = queryData<readType>(nameAtt,reader,sizex,sizey,ischar);
   if(reader == NULL)return ERROR;
   value = *reader;

   if(ischar)deleteifchar(reader);
   delete [] reader;

   return out;
}

template<typename outType, typename readType>
inline
status HDF5Reader::queryVectorData(const std::string &nameAtt, outType &value) const
{
   readType *reader(NULL);
   unsigned int sizex(0),sizey(0);
   bool ischar(false);
   status out = queryData<readType>(nameAtt,reader,sizex,sizey,ischar);
   if(reader == NULL)return ERROR;
   value.resize(sizex);
   for(unsigned int i = 0; i < sizex; i++)
   {
      value[i] = *(reader + i);
      if(ischar)deleteifchar(reader+i);
   }

   delete [] reader;

   return out;
}

template<typename outType, typename readType>
inline
status HDF5Reader::queryMatrixData(const std::string &nameAtt, outType &value) const
{
   readType *reader(NULL);
   unsigned int sizex(0),sizey(0);
   bool ischar(false);
   status out = queryData<readType>(nameAtt,reader,sizex,sizey,ischar);
   if(reader == NULL)return ERROR;
   value.resize(sizex);
   for(unsigned int i = 0; i < sizex; i++)
   {
     value[i].resize(sizey);
     for(unsigned int j = 0; j < sizey; j++)
     {
      value[i][j] = *(reader + j + sizex*i);
      if(ischar)deleteifchar(reader + j + sizex*i);
     }
   }

   delete [] reader;

   return out;
}

template <typename readType>
inline
status HDF5Reader::queryAttribute(const std::string &nameAtt, readType *&reader, 
                                              unsigned int &sizex, unsigned int &sizey, 
                                              bool &ischar, const hid_t &where) const
{
   if(!H5Aexists(where,nameAtt.c_str()))return NOT_FOUND;

   hid_t att_id       = H5Aopen(where,nameAtt.c_str(),H5P_DEFAULT);
   hid_t att_space_id = H5Aget_space(att_id); // set space to obtain sizes
   hid_t att_type_id  = H5Aget_type(att_id);

   int nrank = H5Sget_simple_extent_ndims(att_space_id); // number of rows (values per parameters)
   if(nrank > 2)return ERROR;
   hsize_t dims[nrank]; // nrankAtt, dimension of tensor
   herr_t hstatus = H5Sget_simple_extent_dims(att_space_id,dims,NULL); // fills with each tensor dimension's dimension
   if(hstatus < 0)return ERROR;
   unsigned int fullSize(dims[0]),width(0);
   if(nrank == 2)
   {
     fullSize *= dims[1];
     width = dims[1];
   }
   sizex = dims[0];
   sizey = width;
   reader = new readType[fullSize];

//dims[1]*y + x
//in case of string
   hid_t string_type = H5Tcopy(H5T_C_S1);
   if(H5Tget_class(att_type_id) == H5T_STRING)
   {
     ischar = true;
     if(H5Tis_variable_str(att_type_id))
     {
       string_type = H5Tget_native_type(att_type_id, H5T_DIR_ASCEND);
     }else
     {
       hstatus = H5Tset_size(string_type,H5Tget_size(att_type_id));
       if(hstatus < 0)return ERROR;
     }
   }

   hstatus = H5Aread(att_id,att_type_id,reader);
   if(hstatus < 0)return ERROR;

   if(H5Tget_class(att_type_id) == H5T_STRING)H5Tclose(string_type);

//close
   H5Tclose(att_type_id);
   H5Aclose(att_id);
   H5Sclose(att_space_id);

   return SUCCESS;

}

template<typename outType, typename readType>
inline
status HDF5Reader::queryScalarAttribute(const std::string &nameAtt, outType &value, hid_t where) const
{
   readType *reader(NULL);
   unsigned int sizex(0),sizey(0);
   bool ischar(false);
   if(where == 0)where = me;
   status out = queryAttribute<readType>(nameAtt,reader,sizex,sizey,ischar,where);
   if(reader == NULL)return ERROR;
   value = *reader;

   if(ischar)deleteifchar(reader);
   delete [] reader;

   return out;
}

template<typename outType, typename readType>
inline
status HDF5Reader::queryVectorAttribute(const std::string &nameAtt, outType &value, hid_t where) const
{
   readType *reader(NULL);
   unsigned int sizex(0),sizey(0);
   bool ischar(false);
   if(where == 0)where = me;
   status out = queryAttribute<readType>(nameAtt,reader,sizex,sizey,ischar,where);
   if(reader == NULL)return ERROR;
   value.resize(sizex);
   for(unsigned int i = 0; i < sizex; i++)
   {
      value[i] = *(reader + i);
      if(ischar)deleteifchar(reader+i);
   }

   delete [] reader;

   return out;
}

template<typename outType, typename readType>
inline
status HDF5Reader::queryMatrixAttribute(const std::string &nameAtt, outType &value, hid_t where) const
{
   readType *reader(NULL);
   unsigned int sizex(0),sizey(0);
   bool ischar(false);
   if(where == 0)where = me;
   status out = queryAttribute<readType>(nameAtt,reader,sizex,sizey,ischar,where);
   if(reader == NULL)return ERROR;
   value.resize(sizex);
   for(unsigned int i = 0; i < sizex; i++)
   {
     value[i].resize(sizey);
     for(unsigned int j = 0; j < sizey; j++)
     {
      value[i][j] = *(reader + sizex*i + j);
      if(ischar)deleteifchar(reader + j + sizex*i);
     }
   }

   delete [] reader;

   return out;
}
}//namespace root
}//namespace Antioch
#endif
