//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-
#ifndef MSDATABASEREADER_H
#define MSDATABASEREADER_H

#include <string>
#include <vector>

namespace Antioch
{

namespace root
{

typedef std::vector<double>       VectorDouble;
typedef std::vector<int>          VectorInt;
typedef std::vector<bool>         VectorBool;
typedef std::vector<std::string>  VectorString;
typedef std::vector<VectorDouble> MatrixDouble;
typedef std::vector<VectorInt>    MatrixInt;
typedef std::vector<VectorBool>   MatrixBool;
typedef std::vector<VectorString> MatrixString;
enum status { SUCCESS , WRONG_TYPE , NOT_FOUND , ERROR};

template <typename typeDB, typename idType = std::string>
class DBNodeReader
{
    public:

    DBNodeReader(){};
    DBNodeReader(const typeDB &inter):dbReader(inter){}
    ~DBNodeReader(){dbReader.close();}

    bool load(const std::string & filename) {return dbReader.load(filename);}
    void unload() {dbReader.close();dbReader.unload();}

    void setObj(const typeDB &inter){dbReader = inter;} 

    const idType getId()  const     {return dbReader.getId();}

//this attributes
    status queryDoubleAttribute(const std::string &nameAtt, double &att)     const {return dbReader.queryDoubleAttribute(nameAtt,att);}
    status queryIntAttribute   (const std::string &nameAtt, int &att)        const {return dbReader.queryIntAttribute(nameAtt,att);}
    status queryBoolAttribute  (const std::string &nameAtt, bool &att)       const {return dbReader.queryBoolAttribute(nameAtt,att);}
    status queryStringAttribute(const std::string &nameAtt, std::string &att)const {return dbReader.queryStringAttribute(nameAtt,att);}

    status queryVectorDoubleAttribute(const std::string &nameAtt, VectorDouble &att)const {return dbReader.queryVectorDoubleAttribute(nameAtt,att);}
    status queryVectorIntAttribute   (const std::string &nameAtt, VectorInt &att)   const {return dbReader.queryVectorIntAttribute(nameAtt,att);}
    status queryVectorBoolAttribute  (const std::string &nameAtt, VectorBool &att)  const {return dbReader.queryVectorBoolAttribute(nameAtt,att);}
    status queryVectorStringAttribute(const std::string &nameAtt, VectorString &att)const {return dbReader.queryVectorStringAttribute(nameAtt,att);}

    status queryMatrixDoubleAttribute(const std::string &nameAtt, MatrixDouble &att)const {return dbReader.queryMatrixDoubleAttribute(nameAtt,att);}
    status queryMatrixIntAttribute   (const std::string &nameAtt, MatrixInt &att)   const {return dbReader.queryMatrixIntAttribute(nameAtt,att);}
    status queryMatrixBoolAttribute  (const std::string &nameAtt, MatrixBool &att)  const {return dbReader.queryMatrixBoolAttribute(nameAtt,att);}
    status queryMatrixStringAttribute(const std::string &nameAtt, MatrixString &att)const {return dbReader.queryMatrixStringAttribute(nameAtt,att);}

//this data
    status queryDoubleData(const std::string &nameData, double &data)     const {return dbReader.queryDoubleData(nameData,data);}
    status queryIntData   (const std::string &nameData, int &data)        const {return dbReader.queryIntData(nameData,data);}
    status queryBoolData  (const std::string &nameData, bool &data)       const {return dbReader.queryBoolData(nameData,data);}
    status queryStringData(const std::string &nameData, std::string &data)const {return dbReader.queryStringData(nameData,data);}

    status queryVectorDoubleData(const std::string &nameData, VectorDouble &data)const {return dbReader.queryVectorDoubleData(nameData,data);}
    status queryVectorIntData   (const std::string &nameData, VectorInt &data)   const {return dbReader.queryVectorIntData(nameData,data);}
    status queryVectorBoolData  (const std::string &nameData, VectorBool &data)  const {return dbReader.queryVectorBoolData(nameData,data);}
    status queryVectorStringData(const std::string &nameData, VectorString &data)const {return dbReader.queryVectorStringData(nameData,data);}

    status queryMatrixDoubleData(const std::string &nameData, MatrixDouble &data)const {return dbReader.queryMatrixDoubleData(nameData,data);}
    status queryMatrixIntData   (const std::string &nameData, MatrixInt &data)   const {return dbReader.queryMatrixIntData(nameData,data);}
    status queryMatrixBoolData  (const std::string &nameData, MatrixBool &data)  const {return dbReader.queryMatrixBoolData(nameData,data);}
    status queryMatrixStringData(const std::string &nameData, MatrixString &data)const {return dbReader.queryMatrixStringData(nameData,data);}

//this data's attribute
    status queryDoubleDataAttribute(const std::string &nameData, const std::string &nameAtt, double &att)      const
                                                                                {return dbReader.queryDoubleDataAttribute(nameData,nameAtt,att);}
    status queryIntDataAttribute   (const std::string &nameData, const std::string &nameAtt, int &att)         const
                                                                                {return dbReader.queryIntDataAttribute(nameData,nameAtt,att);}
    status queryBoolDataAttribute  (const std::string &nameData, const std::string &nameAtt, bool &att)        const
                                                                                {return dbReader.queryBoolDataAttribute(nameData,nameAtt,att);}
    status queryStringDataAttribute(const std::string &nameData, const std::string &nameAtt, std::string &att) const 
                                                                                {return dbReader.queryStringDataAttribute(nameData,nameAtt,att);}

    status queryVectorDoubleDataAttribute(const std::string &nameData, const std::string &nameAtt, VectorDouble &att) const 
                                                                                {return dbReader.queryVectorDoubleDataAttribute(nameData,nameAtt,att);}
    status queryVectorIntDataAttribute   (const std::string &nameData, const std::string &nameAtt, VectorInt &att)    const 
                                                                                {return dbReader.queryVectorIntDataAttribute(nameData,nameAtt,att);}
    status queryVectorBoolDataAttribute  (const std::string &nameData, const std::string &nameAtt, VectorBool &att)   const 
                                                                                {return dbReader.queryVectorBoolDataAttribute(nameData,nameAtt,att);}
    status queryVectorStringDataAttribute(const std::string &nameData, const std::string &nameAtt, VectorString &att) const 
                                                                                {return dbReader.queryVectorStringDataAttribute(nameData,nameAtt,att);}

    status queryMatrixDoubleDataAttribute(const std::string &nameData, const std::string &nameAtt, MatrixDouble &att) const 
                                                                                {return dbReader.queryMatrixDoubleDataAttribute(nameData,nameAtt,att);}
    status queryMatrixIntDataAttribute   (const std::string &nameData, const std::string &nameAtt, MatrixInt &att)    const 
                                                                                {return dbReader.queryMatrixIntDataAttribute(nameData,nameAtt,att);}
    status queryMatrixBoolDataAttribute  (const std::string &nameData, const std::string &nameAtt, MatrixBool &att)   const 
                                                                                {return dbReader.queryMatrixBoolDataAttribute(nameData,nameAtt,att);}
    status queryMatrixStringDataAttribute(const std::string &nameData, const std::string &nameAtt, MatrixString &att) const 
                                                                                {return dbReader.queryMatrixStringDataAttribute(nameData,nameAtt,att);}

    DBNodeReader<typeDB,idType>* firstChild() const 
                              {
                                DBNodeReader<typeDB,idType> *child = new DBNodeReader<typeDB,idType>(dbReader.firstChild());
                                return child;
                              }
    DBNodeReader<typeDB,idType>* nextSibling()const 
                              {
                                DBNodeReader<typeDB,idType> *sibling = new DBNodeReader<typeDB,idType>(dbReader.nextSibling());
                                return sibling;
                              }

    typeDB dbReader;
 };


template <typename typeDB, typename idType = std::string>
class DBNodeWriter
{
    public:

    DBNodeWriter(){};
    DBNodeWriter(const typeDB &inter):dbWriter(inter){}
    ~DBNodeWriter(){};

    void setObj(const typeDB &inter){dbWriter = inter;} 

    bool load(const std::string & filename) {return dbWriter.load(filename);}

    const idType getId()  const     {return dbWriter.getId();}

//this attributes
    status writeDoubleAttribute(const std::string &nameAtt, const double &att)     {return dbWriter.writeDoubleAttribute(nameAtt,att);}
    status writeIntAttribute   (const std::string &nameAtt, const int &att)        {return dbWriter.writeIntAttribute(nameAtt,att);}
    status writeBoolAttribute  (const std::string &nameAtt, const bool &att)       {return dbWriter.writeBoolAttribute(nameAtt,att);}
    status writeStringAttribute(const std::string &nameAtt, const std::string &att){return dbWriter.writeStringAttribute(nameAtt,att);}

    status writeVectorDoubleAttribute(const std::string &nameAtt, const VectorDouble &att){return dbWriter.writeVectorDoubleAttribute(nameAtt,att);}
    status writeVectorIntAttribute   (const std::string &nameAtt, const VectorInt &att)   {return dbWriter.writeVectorIntAttribute(nameAtt,att);}
    status writeVectorBoolAttribute  (const std::string &nameAtt, const VectorBool &att)  {return dbWriter.writeVectorBoolAttribute(nameAtt,att);}
    status writeVectorStringAttribute(const std::string &nameAtt, const VectorString &att){return dbWriter.writeVectorStringAttribute(nameAtt,att);}

    status writeMatrixDoubleAttribute(const std::string &nameAtt, const MatrixDouble &att){return dbWriter.writeMatrixDoubleAttribute(nameAtt,att);}
    status writeMatrixIntAttribute   (const std::string &nameAtt, const MatrixInt &att)   {return dbWriter.writeMatrixIntAttribute(nameAtt,att);}
    status writeMatrixBoolAttribute  (const std::string &nameAtt, const MatrixBool &att)  {return dbWriter.writeMatrixBoolAttribute(nameAtt,att);}
    status writeMatrixStringAttribute(const std::string &nameAtt, const MatrixString &att){return dbWriter.writeMatrixStringAttribute(nameAtt,att);}

//this data
    status writeDoubleData(const std::string &nameData, const double &data)     {return dbWriter.writeDoubleData(nameData,data);}
    status writeIntData   (const std::string &nameData, const int &data)        {return dbWriter.writeIntData(nameData,data);}
    status writeBoolData  (const std::string &nameData, const bool &data)       {return dbWriter.writeBoolData(nameData,data);}
    status writeStringData(const std::string &nameData, const std::string &data){return dbWriter.writeStringData(nameData,data);}

    status writeVectorDoubleData(const std::string &nameData, const VectorDouble &data){return dbWriter.writeVectorDoubleData(nameData,data);}
    status writeVectorIntData   (const std::string &nameData, const VectorInt &data)   {return dbWriter.writeVectorIntData(nameData,data);}
    status writeVectorBoolData  (const std::string &nameData, const VectorBool &data)  {return dbWriter.writeVectorBoolData(nameData,data);}
    status writeVectorStringData(const std::string &nameData, const VectorString &data){return dbWriter.writeVectorStringData(nameData,data);}

    status writeMatrixDoubleData(const std::string &nameData, const MatrixDouble &data){return dbWriter.writeMatrixDoubleData(nameData,data);}
    status writeMatrixIntData   (const std::string &nameData, const MatrixInt &data)   {return dbWriter.writeMatrixIntData(nameData,data);}
    status writeMatrixBoolData  (const std::string &nameData, const MatrixBool &data)  {return dbWriter.writeMatrixBoolData(nameData,data);}
    status writeMatrixStringData(const std::string &nameData, const MatrixString &data){return dbWriter.writeMatrixStringData(nameData,data);}

//this data's attribute
    status writeDoubleDataAttribute(const std::string &nameData, const std::string &nameAtt, const double &att)
                                                                                {return dbWriter.writeDoubleDataAttribute(nameData,nameAtt,att);}
    status writeIntDataAttribute   (const std::string &nameData, const std::string &nameAtt, const int &att)
                                                                                {return dbWriter.writeIntDataAttribute(nameData,nameAtt,att);}
    status writeBoolDataAttribute  (const std::string &nameData, const std::string &nameAtt, const bool &att)
                                                                                {return dbWriter.writeBoolDataAttribute(nameData,nameAtt,att);}
    status writeStringDataAttribute(const std::string &nameData, const std::string &nameAtt, const std::string &att)
                                                                                {return dbWriter.writeStringDataAttribute(nameData,nameAtt,att);}

    status writeVectorDoubleDataAttribute(const std::string &nameData, const std::string &nameAtt, const VectorDouble &att)
                                                                                {return dbWriter.writeVectorDoubleDataAttribute(nameData,nameAtt,att);}
    status writeVectorIntDataAttribute   (const std::string &nameData, const std::string &nameAtt, const VectorInt &att)
                                                                                {return dbWriter.writeVectorIntDataAttribute(nameData,nameAtt,att);}
    status writeVectorBoolDataAttribute  (const std::string &nameData, const std::string &nameAtt, const VectorBool &att)
                                                                                {return dbWriter.writeVectorBoolDataAttribute(nameData,nameAtt,att);}
    status writeVectorStringDataAttribute(const std::string &nameData, const std::string &nameAtt, const VectorString &att)
                                                                                {return dbWriter.writeVectorStringDataAttribute(nameData,nameAtt,att);}

    status writeMatrixDoubleDataAttribute(const std::string &nameData, const std::string &nameAtt, const MatrixDouble &att)
                                                                                {return dbWriter.writeMatrixDoubleDataAttribute(nameData,nameAtt,att);}
    status writeMatrixIntDataAttribute   (const std::string &nameData, const std::string &nameAtt, const MatrixInt &att)
                                                                                {return dbWriter.writeMatrixIntDataAttribute(nameData,nameAtt,att);}
    status writeMatrixBoolDataAttribute  (const std::string &nameData, const std::string &nameAtt, const MatrixBool &att)
                                                                                {return dbWriter.writeMatrixBoolDataAttribute(nameData,nameAtt,att);}
    status writeMatrixStringDataAttribute(const std::string &nameData, const std::string &nameAtt, const MatrixString &att)
                                                                                {return dbWriter.writeMatrixStringDataAttribute(nameData,nameAtt,att);}

    DBNodeWriter<typeDB, idType> * addChild(const idType &value) {DBNodeWriter<typeDB,idType> *child = new DBNodeWriter(dbWriter.addChild(value)); 
                                                                       return child;}
    typeDB dbWriter;
 };

template <typename typeLoad>
class DBTreeLoader
{
    public:

    DBTreeLoader():root(NULL){}
    DBTreeLoader(const std::string &filename):root(NULL),file(filename){if(!load())delete root;}
    ~DBTreeLoader(){unload();}

    bool load()  {
                   if(root != NULL)delete root;
                   root = new typeLoad;
                   return root->load(file);
                 }

    void unload(){if(root != NULL){root->unload();delete root;}}

    bool load(const std::string &filename){
                                           unload();
                                           file = filename;
                                           return load();
                                          }

    typeLoad * getRootNode()const {return root;}

    private:

    typeLoad * root;
    std::string file;
};

template <typename typeDB, typename idType = std::string>
struct DBTreeWriter
{
  typedef DBTreeLoader<DBNodeWriter<typeDB, idType> > write;
};

template <typename typeDB, typename idType = std::string>
struct DBTreeReader
{
  typedef DBTreeLoader<DBNodeReader<typeDB, idType> > read;
};

}//namespace root
}//namespace Antioch
#endif // MSDATABASEREADER_H
