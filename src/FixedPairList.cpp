#include <sstream>
#include "python.hpp"
#include "FixedPairList.hpp"
//#include <vector>
//#include <utility>
//#include <algorithm>
#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "Buffer.hpp"

#include "esutil/Error.hpp"

using namespace std;

namespace espresso {

  /*
  FixedPairList::FixedPairList(shared_ptr< storage::Storage > _storage)
  : FixedListComm (_storage){}

  FixedPairList::~FixedPairList() {
    //std::cout << "~fixedpairlist" << std::endl;
    //FixedListComm::~FixedListComm();
  }*/


  LOG4ESPP_LOGGER(FixedPairList::theLogger, "FixedPairList");


  FixedPairList::FixedPairList(shared_ptr< storage::Storage > _storage)
    : storage(_storage), globalPairs()
  {
    LOG4ESPP_INFO(theLogger, "construct FixedPairList");

    con1 = storage->beforeSendParticles.connect
      (boost::bind(&FixedPairList::beforeSendParticles, this, _1, _2));
    con2 = storage->afterRecvParticles.connect
      (boost::bind(&FixedPairList::afterRecvParticles, this, _1, _2));
    con3 = storage->onParticlesChanged.connect
      (boost::bind(&FixedPairList::onParticlesChanged, this));
  }

  FixedPairList::~FixedPairList() {

    LOG4ESPP_INFO(theLogger, "~FixedPairList");

    con1.disconnect();
    con2.disconnect();
    con3.disconnect();
  }


  /*
  bool FixedPairList::add(longint pid1, longint pid2) {
    std::vector<longint> tmp;
    tmp.push_back(pid2);
    tmp.push_back(pid1); // this is used as key

    return FixedListComm::add(tmp);
  }*/


  bool FixedPairList::
  add(longint pid1, longint pid2) {
    bool returnVal = true;
    if (pid1 > pid2)
      std::swap(pid1, pid2);

    System& system = storage->getSystemRef();
    esutil::Error err(system.comm);
    
    // ADD THE LOCAL PAIR
    Particle *p1 = storage->lookupRealParticle(pid1);
    Particle *p2 = storage->lookupLocalParticle(pid2);
    
    if (!p1){
      // Particle does not exist here, return false
      returnVal=false;
    }
    else{
      if (!p2) {
        std::stringstream msg;
        msg << "Adding error. Fixed Pair List particle p2 " << pid2 << " does not exists here and cannot be added";
        err.setException( msg.str() );
        //std::runtime_error(err.str());
      }
    }
    err.checkException();
    
    if(returnVal){
      // add the pair locally
      this->add(p1, p2);

      // add the particle pair to the globalPairs list
      globalPairs.insert(std::make_pair(pid1, pid2));
    }
    LOG4ESPP_INFO(theLogger, "added fixed pair to global pair list");
    return returnVal;
  }

  python::list FixedPairList::getBonds()
  {
	python::tuple bond;
	python::list bonds;
	for (GlobalPairs::const_iterator it=globalPairs.begin(); it != globalPairs.end(); it++) {
      bond = python::make_tuple(it->first, it->second);
      bonds.append(bond);
    }

	return bonds;
  }

  void FixedPairList::
  beforeSendParticles(ParticleList& pl, 
		      OutBuffer& buf) {
    std::vector< longint > toSend;
    // loop over the particle list
    for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
      longint pid = pit->id();
      
      LOG4ESPP_DEBUG(theLogger, "send particle with pid " << pid << ", find pairs");

      // find all pairs that involve this particle
      
      int n = globalPairs.count(pid);

      if (n > 0) {
        std::pair<GlobalPairs::const_iterator,
          GlobalPairs::const_iterator> equalRange
          = globalPairs.equal_range(pid);

        // first write the pid of the first particle
        // then the number of partners
        // and then the pids of the partners
        toSend.reserve(toSend.size()+n+1);
        toSend.push_back(pid);
        toSend.push_back(n);
        for (GlobalPairs::const_iterator it = equalRange.first;
             it != equalRange.second; ++it) {
          toSend.push_back(it->second);
          LOG4ESPP_DEBUG(theLogger, "send global bond: pid "
                       << pid << " and partner " << it->second);
        }

        // delete all of these pairs from the global list
        //globalPairs.erase(equalRange.first->first, equalRange.second->first);
        globalPairs.erase(pid);
        // std::cout << "erasing particle " << pid << " from here" << std::endl;
      }
    }
    // send the list
    buf.write(toSend);
    LOG4ESPP_INFO(theLogger, "prepared fixed pair list before send particles");
  }

  void FixedPairList::
  afterRecvParticles(ParticleList &pl, 
		     InBuffer& buf) {
    std::vector< longint > received;
    int n;
    longint pid1, pid2;
    GlobalPairs::iterator it = globalPairs.begin();
    // receive the bond list
    buf.read(received);
    int size = received.size(); int i = 0;
    while (i < size) {
      // unpack the list
      pid1 = received[i++];
      n = received[i++];
      LOG4ESPP_DEBUG(theLogger, "recv particle " << pid1 << 
                                ", has " << n << " global pairs");
      for (; n > 0; --n) {
	pid2 = received[i++];
	// add the bond to the global list
        LOG4ESPP_DEBUG(theLogger, "received pair " << pid1 << " , " << pid2);
	it = globalPairs.insert(it, std::make_pair(pid1, pid2));
      }
    }
    if (i != size) {
      LOG4ESPP_ERROR(theLogger, 
        "ATTETNTION:  recv particles might have read garbage\n");
    }
    LOG4ESPP_INFO(theLogger, "received fixed pair list after receive particles");
  }

  void FixedPairList::
  onParticlesChanged() {
    LOG4ESPP_INFO(theLogger, "rebuild local bond list from global\n");

    System& system = storage->getSystemRef();
    esutil::Error err(system.comm);
    
    this->clear();
    longint lastpid1 = -1;
    Particle *p1;
    Particle *p2;
    for (GlobalPairs::const_iterator it = globalPairs.begin();
	 it != globalPairs.end(); ++it) {
      if (it->first != lastpid1) {
	    p1 = storage->lookupRealParticle(it->first);
        if (p1 == NULL) {
          std::stringstream msg;
          msg << "onParticlesChanged error. Fixed Pair List particle p1 " << it->first << " does not exists here";
          err.setException( msg.str() );
          //std::runtime_error(err.str());
        }
	    lastpid1 = it->first;
      }
      p2 = storage->lookupLocalParticle(it->second);
      if (p2 == NULL) {
          std::stringstream msg;
          msg << "onParticlesChanged error. Fixed Pair List particle p2 " << it->second << " does not exists here";
          //std::runtime_error(err.str());
          err.setException( msg.str() );
      }
      this->add(p1, p2);
    }
    err.checkException();
    
    LOG4ESPP_INFO(theLogger, "regenerated local fixed pair list from global list");
  }

  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/

  void FixedPairList::registerPython() {

    using namespace espresso::python;

    bool (FixedPairList::*pyAdd)(longint pid1, longint pid2)
      = &FixedPairList::add;
    //bool (FixedPairList::*pyAdd)(pvec pids) = &FixedPairList::add;

    class_<FixedPairList, shared_ptr<FixedPairList> >
      ("FixedPairList", init <shared_ptr<storage::Storage> >())
      .def("add", pyAdd)
      .def("size", &FixedPairList::size)
      .def("getBonds",  &FixedPairList::getBonds)
      ;
  }
}
