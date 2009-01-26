#include "PHOTOSintf.h"
#include <vector>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <iostream>

class LHEF {
public:

  struct Event {
    int n_event;
    int nup, idprup;
    double xwgtup, scalup, aqedup, aqcdup;
    struct Particle {
      int pdgid;
      int status;
      int mother1,mother2;
      int color1,color2;
      double px,py,pz,E,m;
      double lifetime,spin;
    };
    std::vector<Particle> particles;
  };

  bool loadNextEvent();
  void storeEvent();

  Event event;

  void loadFromFile(const std::string& fname);
  void storeToFile(const std::string& fname, const std::string& pretext);

  const std::string& getPretext() const { return m_pretext; }

  LHEF() : m_file(0) {}

  ~LHEF() {     
    if (m_file!=0) {
      if (m_writing) fprintf(m_file,"</LesHouchesEvents>\n");
      fclose(m_file);
    }
    m_file=0;
  }

private:
  std::string m_pretext;
  FILE* m_file;
  bool m_writing;
};

void event2hepevt(LHEF::Event& evt);
void hepevt2event(LHEF::Event& evt,bool);

int main(int argc, char* argv[]) {

    if (argc<3 || argc>4)
    {
        std::cout << "usage: photosLHEF  file_in  file_out  [1/0 i.e. keep/nokeep]" << std::endl;
        exit (0);
    }

  phoini_();

  int mod=3;
  phcork_(mod);

  LHEF lhefIn, lhefOut;  
  
  // load event from LHEF
  lhefIn.loadFromFile(argv[1]);
  lhefOut.storeToFile(argv[2],lhefIn.getPretext());

  bool keep(false);
  if(argc==4) keep=atoi(argv[3]); 
  
  while (lhefIn.loadNextEvent()) {
    event2hepevt(lhefIn.event);
    
    bool showMe=false;

    if (lhefIn.event.n_event<3) showMe=true;
    if (lhefIn.event.n_event==439) showMe=true;
    
    if (showMe)
      phodmp_();
    int which=1;
    photos_(which);
    if (showMe)
      phodmp_();

    lhefOut.event=lhefIn.event;
    hepevt2event(lhefOut.event,keep);
    lhefOut.storeEvent();
  }

  return 0;
}


void event2hepevt(LHEF::Event& evt) {
  hepevt_.d_h_nevhep=evt.n_event;
  hepevt_.d_h_nhep=evt.particles.size();
  int j=0;
  for (std::vector<LHEF::Event::Particle>::const_iterator i=evt.particles.begin();
       i!=evt.particles.end(); ++i) {
    hepevt_.d_h_isthep[j]=std::max(0,i->status);
    hepevt_.d_h_idhep[j]=i->pdgid;
    hepevt_.d_h_jmohep[j][0]=i->mother1;
    hepevt_.d_h_jmohep[j][1]=i->mother2;

    hepevt_.d_h_phep[j][0]=i->px;
    hepevt_.d_h_phep[j][1]=i->py;
    hepevt_.d_h_phep[j][2]=i->pz;
    hepevt_.d_h_phep[j][3]=i->E;
    hepevt_.d_h_phep[j][4]=i->m;

    if (i->m==0 && abs(i->pdgid)==11)  hepevt_.d_h_phep[j][4]=0.0005;
    if (i->m==0 && abs(i->pdgid)==2)  hepevt_.d_h_phep[j][4]=0.010;
    if (i->m==0 && abs(i->pdgid)==1)  hepevt_.d_h_phep[j][4]=0.010;
    
    hepevt_.d_h_vhep[j][0]=0;
    hepevt_.d_h_vhep[j][1]=0;
    hepevt_.d_h_vhep[j][2]=0;
    hepevt_.d_h_vhep[j][3]=0;
    ++j;
  }

  for (int j=0; j<int(evt.particles.size()); j++) {
    // daughters 
    int d1=0,d2=0;
    for (int jj=j+1; jj<int(evt.particles.size()); jj++) {
      if (hepevt_.d_h_jmohep[jj][0]==j+1 || (hepevt_.d_h_jmohep[jj][1]==j+1)) {
	if (d1==0) d1=jj+1;
	else d2=jj+1;
      }
    }
    if (d2==0) d2=d1;
    hepevt_.d_h_jdahep[j][0]=d1;
    hepevt_.d_h_jdahep[j][1]=d2;
  }
}

void hepevt2event(LHEF::Event& evt, bool keep) {
  // just do the particles beyond the end...
  // int had=evt.particles.size();

  // clear and reinsert all particles in event
  // (not just those added by photos) to allow flagging of 'too much brem'-like events
  evt.particles.clear();

  for (int j=0; j<hepevt_.d_h_nhep; ++j) {
    LHEF::Event::Particle p;
    p.pdgid=hepevt_.d_h_idhep[j];
    p.mother1=hepevt_.d_h_jmohep[j][0];
    p.mother2=hepevt_.d_h_jmohep[j][1];
    
    p.status=hepevt_.d_h_isthep[j];

    // if there's a particle with 'too much brem'-like status (9999)
    // drop or keep event according to user's flag. Drop is default.
    if(p.status==9999 && (!keep))
    {
        evt.particles.clear();
        return;
    }

    p.color1=0;
    p.color2=0;

    p.px=hepevt_.d_h_phep[j][0];
    p.py=hepevt_.d_h_phep[j][1];
    p.pz=hepevt_.d_h_phep[j][2];
    p.E=hepevt_.d_h_phep[j][3]; 
    p.m=hepevt_.d_h_phep[j][4]; 

    p.lifetime=0;
    p.spin=0;

    evt.particles.push_back(p);
  }
}
 
bool LHEF::loadNextEvent() {
  char buffer[1024];
  std::vector<std::string> lines;
  int mode=0;
  while (!feof(m_file)) {
    buffer[0]=0;
    fgets(buffer,1000,m_file);
    if (strstr(buffer,"<event>")) mode=1;
    else if (strstr(buffer,"</event>")) {
      mode=2;
      break;
    } else if (mode==1) lines.push_back(buffer);
  }
  if (mode!=2 || lines.empty()) return false;
  event.n_event++;
  event.particles.clear();
  int got;
  got=sscanf(lines[0].c_str()," %d %d %lf %lf %lf %lf ",&event.nup,&event.idprup,
	     &event.xwgtup, &event.scalup, &event.aqedup, &event.aqcdup);
  if (got!=6) {
    fprintf(stderr,"Parse error: Looking for 6, saw %d\n",got);
    return false;
  }
  if (int(lines.size())!=event.nup+1) {
    fprintf(stderr,"Parse error: Expected %d lines, saw %d\n",event.nup+1,lines.size());
    return false;
  }
  for (int i=1; i<=event.nup; i++) {
    Event::Particle p;
    got=sscanf(lines[i].c_str()," %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf ",&p.pdgid,&p.status,&p.mother1,&p.mother2,&p.color1,&p.color2,&p.px,&p.py,&p.pz,&p.E,&p.m,&p.lifetime,&p.spin);
    if (got!=6+7) {
      fprintf(stderr,"Parse error: Line %d Looking for 6+7, saw %d\n",i,got);
      return false;
    }
    event.particles.push_back(p);
  }
  return true;
  
}

void LHEF::storeEvent() {
  fprintf(m_file,"<event>\n");
  fprintf(m_file ," %d %d %le %le %le %le\n",event.particles.size(),event.idprup,
	  event.xwgtup,event.scalup,event.aqedup,event.aqcdup);

  for (std::vector<Event::Particle>::const_iterator i=event.particles.begin();
       i!=event.particles.end(); ++i) {
    
    fprintf(m_file," %3d %2d %2d %2d %3d %3d %+le %+le %+le %+le %+le %.1le %.2le\n",i->pdgid,i->status,i->mother1,i->mother2,i->color1,i->color2,i->px,i->py,i->pz,i->E,i->m,i->lifetime,i->spin);
  }
  fprintf(m_file,"</event>\n");
}

void LHEF::loadFromFile(const std::string& fname) {
  m_writing=false;
  m_file=fopen(fname.c_str(),"r");
  event.n_event=0;

  char buffer[1024];
  bool done=false;
  while (!done && !feof(m_file)) {
    buffer[0]=0;
    fgets(buffer,1000,m_file);
    if (strstr(buffer,"</init>")) done=true;
    m_pretext+=buffer;
  }

}

void LHEF::storeToFile(const std::string& fname, const std::string& pretext) {
  m_writing=true;
  m_file=fopen(fname.c_str(),"w");
  event.n_event=0;
  fputs(pretext.c_str(),m_file);
}
