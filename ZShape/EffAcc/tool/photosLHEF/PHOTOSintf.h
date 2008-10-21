extern "C"
{
  extern struct {
    int d_h_nevhep;
    int d_h_nhep;
    int d_h_isthep[4000];
    int d_h_idhep[4000];
    int d_h_jmohep[4000][2];
    int d_h_jdahep[4000][2];
    double d_h_phep[4000][5];
    double d_h_vhep[4000][4];
  } hepevt_;

  extern struct {
    double alpha;
    double xphcut;
  } phocop_;
}

extern "C" {
  void phoini_(); // initialize
  void photos_get_(); // copy into local data
  void photos_make_(int& id); // 
  void photos_set_(); // copy from local data
  void photos_(int& id); // 
  void phodmp_(); // dump
  void phcork_(int &mod); //
}
