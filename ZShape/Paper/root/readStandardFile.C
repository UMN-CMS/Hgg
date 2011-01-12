TH1* readStandardFile(const char* name, const char* fname) {
  TH1* dist=new TH1F(name,name,100,-5.0,5.0);
  dist->SetDirectory(0);
  dist->Sumw2();

  FILE* af=fopen(fname,"rt");
  char buffer[1024];

  int irow;
  float flow,fhi,val,err;
  while (!feof(af)) {
    buffer[0]=0;
    fgets(buffer,1023,af);
    if (strchr(buffer,'#')!=0) *(strchr(buffer,'#'))=0;
    
    int matched=sscanf(buffer,"%d %f %f %f %f",&irow,&flow,&fhi,&val,&err);

    dist->SetBinContent(irow,val);
    dist->SetBinError(irow,err);
  }
  fclose(af);

  return dist;
}
