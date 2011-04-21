
// structure for holding the data for one timestamp, i.e. one row from the data
// file
struct DataRow {

};

// class for storing, reading and manipulating the data
// generally speaking the amounts of data are huge, it should be taken into
// account here
class Data {
    vector<DataRow> data; // vector of data rows, one row is one timestamp
  public:
    Data(); // constructor
};

