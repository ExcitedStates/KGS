/*
    KGSX: Biomolecular Kino-geometric Sampling and Fitting of Experimental Data
    Yao et al, Proteins. 2012 Jan;80(1):25-43
    e-mail: latombe@cs.stanford.edu, vdbedem@slac.stanford.edu, julie.bernauer@inria.fr

        Copyright (C) 2011-2013 Stanford University

        Permission is hereby granted, free of charge, to any person obtaining a copy of
        this software and associated documentation files (the "Software"), to deal in
        the Software without restriction, including without limitation the rights to
        use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
        of the Software, and to permit persons to whom the Software is furnished to do
        so, subject to the following conditions:

        This entire text, including the above copyright notice and this permission notice
        shall be included in all copies or substantial portions of the Software.

        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
        IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
        FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
        OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
        FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
        IN THE SOFTWARE.


*/

#include <assert.h>
#include <iostream>
#include <list>
#include <cmath>

#include "SugarVertex.h"

using namespace std;

SugarVertex::SugarVertex(): RigidbodyGraphVertex(){
    //log("debugRas")<<"SugarVertex() ["<<this<<"]"<<endl;
    DOF_id = -1;
	Cycle_DOF_id = -1;
	isRibose = true;
}

SugarVertex::SugarVertex(int id, Rigidbody* rb): RigidbodyGraphVertex(id,rb){
    //log("debugRas")<<"SugarVertex("<<id<<", "<<rb<<") ["<<this<<"]"<<endl;
	DOF_id = -1;
	Cycle_DOF_id = -1;
	isRibose = true;
}
SugarVertex::~SugarVertex(){}

/* Return true iff the atom is part of a sugar-ring (proline or DNA/RNA ribose). */
bool SugarVertex::sugarMember(Atom* a){
    if(a->getResidue()->getName()=="PRO"){
        if(  a->getName()=="CA" || a->getName()=="CB" || a->getName()=="CG" || a->getName()=="CD" || a->getName()=="N"  ){
            //Check that a sugar-neighbor exists in this RB. Necessary for poly-pro sequences since two Proline CAs can exist in the same RB
            for(vector<Atom*>::iterator nIt = a->Cov_neighbor_list.begin(); nIt!=a->Cov_neighbor_list.end(); nIt++){
                Atom* n = *nIt;
                if( Rb_ptr->containsAtom(n) && (n->getName()=="CA" || n->getName()=="CB" || n->getName()=="CG" || n->getName()=="CD" || n->getName()=="N") )
                    return true;
            }
        }
        return false;
    }
	return 
		(	a->getResidue()->getName()=="G"   &&	(a->getName()=="C1'" || a->getName()=="C2'" || a->getName()=="C3'" || a->getName()=="C4'" || a->getName()=="O4'")	) ||
		(	a->getResidue()->getName()=="C"   &&	(a->getName()=="C1'" || a->getName()=="C2'" || a->getName()=="C3'" || a->getName()=="C4'" || a->getName()=="O4'")	) ||
		(	a->getResidue()->getName()=="A"   &&	(a->getName()=="C1'" || a->getName()=="C2'" || a->getName()=="C3'" || a->getName()=="C4'" || a->getName()=="O4'")	) ||
		(	a->getResidue()->getName()=="T"   &&	(a->getName()=="C1'" || a->getName()=="C2'" || a->getName()=="C3'" || a->getName()=="C4'" || a->getName()=="O4'")	) ||
		(	a->getResidue()->getName()=="U"   &&	(a->getName()=="C1'" || a->getName()=="C2'" || a->getName()=="C3'" || a->getName()=="C4'" || a->getName()=="O4'")	) ;
}

/* Perform a DFS starting at 'a' visiting all connected atoms in 'rb' stopping at sugar-ring atoms. It is assumed 
 * that exactly one sugar atom will be found which is placed at the beginning of the 'visited'-vector.
 * It is also assumed that only one sugar ring is found in 'rb'.*/
void SugarVertex::collectAttached(Atom* a, Rigidbody* rb, list<Atom*>& visited){
	if(std::find(visited.begin(), visited.end(), a)!=visited.end()) return;//Stop if a is visited
	if(!rb->containsAtom(a)) return;	//Stop if a is not a member of this rb
	if(!visited.empty() && sugarMember(visited.front()) && sugarMember(a)) return;//Stop if a sugar-member is already added

	//If a is part of a sugar-ring push it to the front
	if(sugarMember(a)) 	visited.push_front(a);
	else				visited.push_back(a);
	
	//Visit adjacent atoms
	for(vector<Atom*>::iterator nit = a->Cov_neighbor_list.begin(); nit!=a->Cov_neighbor_list.end(); nit++)
		collectAttached(*nit, rb, visited);	

	
}

/* Assuming that 'ringAttached[0].front()' is a sugar-ring atom, this method fills the rest of the 'ringAttached'-array with 
 * atoms such that the first entries form a cycle around the sugar ring and the remaining are connected to the first. */
void SugarVertex::collectRest(Rigidbody* rb){
	for(int i=1;i<5;i++){
		//Atom* prevRingMember = ringAttached[i-1].front();
		for(vector<Atom*>::iterator nit=ringAttached[i-1].front()->Cov_neighbor_list.begin(); nit!=ringAttached[i-1].front()->Cov_neighbor_list.end(); nit++){
			if( (i==1||ringAttached[i-2].front()!=(*nit)) && sugarMember(*nit) ){
				collectAttached(*nit, rb, ringAttached[i]);
				break;
			}
		}
	}

	//Make sure proline and ribose relates the tau DOF to the backbone torsion.
	if(ringAttached[1].front()->getName()=="O4'" || ringAttached[1].front()->getName()=="CD"){
		list<Atom*> tmp = ringAttached[1];
		ringAttached[1] = ringAttached[4];
		ringAttached[4] = tmp;
		tmp = ringAttached[2];
		ringAttached[2] = ringAttached[3];
		ringAttached[3] = tmp;
	}

}

void SugarVertex::setParent(RigidbodyGraphVertex* v){
    //log("debugRas")<<"SugarVertex::setParent("<<v->Rb_ptr<<") this: "<<Rb_ptr<<endl;
	RigidbodyGraphVertex::setParent(v);

	//Determine entry atom
	Atom* entryAtom = NULL;
	map<unsigned int,Edge*>::iterator edge_itr;
	//for (edge_itr=Parent->Edges.begin(); edge_itr != Parent->Edges.end(); ++edge_itr) {
	for (auto eit=Parent->edges.begin(); eit!=Parent->edges.end(); ++edge_itr) {
		//if(RigidbodyGraphVertex::Rb_ptr->containsAtom( (*edge_itr).second->getBond()->Atom2 )) {
		if(RigidbodyGraphVertex::Rb_ptr->containsAtom( (*eit)->getBond()->Atom2 )) {
			//entryAtom = (*edge_itr).second->getBond()->Atom2;
			entryAtom = (*eit)->getBond()->Atom2;
		}
	}
    //Atom* a = Rb_ptr->getAtom("CA");
    //if(a!=NULL) log("debugRas")<<Rb_ptr<<".sugarMember("<<a<<") = "<<sugarMember(a)<<endl;
    collectAttached(entryAtom, Rb_ptr, ringAttached[0] );
    //log("debugRas")<<"Ringattached[0]: ";
    //for(list<Atom*>::iterator aIt = ringAttached[0].begin(); aIt!=ringAttached[0].end(); aIt++){
    //    log("debugRas")<<*aIt<<" , ";
    //}
    //log("debugRas")<<endl;
	collectRest(Rb_ptr);

	Vector3* pos[] = {
		&ringAttached[0].front()->m_Position,
		&ringAttached[1].front()->m_Position,
		&ringAttached[2].front()->m_Position,
		&ringAttached[3].front()->m_Position,
		&ringAttached[4].front()->m_Position
	};

	angles[0] = VectorAngle(*pos[1]-*pos[0],*pos[4]-*pos[0]);
	angles[1] = VectorAngle(*pos[2]-*pos[1],*pos[0]-*pos[1]);
	angles[2] = VectorAngle(*pos[3]-*pos[2],*pos[1]-*pos[2]);
	angles[3] = VectorAngle(*pos[4]-*pos[3],*pos[2]-*pos[3]);
	angles[4] = VectorAngle(*pos[0]-*pos[4],*pos[3]-*pos[4]);
	dists[0] = (*pos[0]-*pos[1]).length();
	dists[1] = (*pos[1]-*pos[2]).length();
	dists[2] = (*pos[2]-*pos[3]).length();
	dists[3] = (*pos[3]-*pos[4]).length();
	dists[4] = (*pos[4]-*pos[0]).length();

	initUp = getCurrentUP();

	//Estimate amplitude of ribose configuration using binary search 
	//Vector3* tmp = NULL;
	Vector3 tmp[20];
	double lower = 0, upper = 3.141592;
	while(upper-lower>0.0001){
		double tmpAm = (lower+upper)/2;
		if(getSugarConformation(0, tmp, tmpAm/*, false*/)) 	{lower = tmpAm;}
		else											{upper = tmpAm;}
        //log("debugRas")<<"Loop .. tested amplitude: "<<tmpAm<<endl;
		//if(lower==tmpAm){
        //	log("debugRas")<<"Sphere["<<(tmp[4*0]).x<<","<<(tmp[4*0]).y<<","<<(tmp[4*0]).z<<",0.3,0.3,0.3,0.3]#C4"<<endl;
        //	log("debugRas")<<"Sphere["<<(tmp[4*1]).x<<","<<(tmp[4*1]).y<<","<<(tmp[4*1]).z<<",0.3,0.3,0.3,0.3]#C3"<<endl;
        //	log("debugRas")<<"Sphere["<<(tmp[4*2]).x<<","<<(tmp[4*2]).y<<","<<(tmp[4*2]).z<<",0.3,0.3,0.3,0.3]#C2"<<endl;
        //	log("debugRas")<<"Sphere["<<(tmp[4*3]).x<<","<<(tmp[4*3]).y<<","<<(tmp[4*3]).z<<",0.3,0.3,0.3,0.3]#C1"<<endl;
        //	log("debugRas")<<"Sphere["<<(tmp[4*4]).x<<","<<(tmp[4*4]).y<<","<<(tmp[4*4]).z<<",0.3,0.8,0.3,0.3]#O4"<<endl;
		//}
	}
	Amplitude = lower;
    //log("debugRas")<<"Sugar amplitude: "<<(Amplitude*180/CTK_PI)<<endl;


	initTorsion = TorsionalAngle(*pos[4], *pos[0], *pos[1], *pos[2]);
	double diff = Amplitude-fabs(initTorsion); 
	if(diff<0.0001 && diff>-0.0001) initTorsion+=Signum(initTorsion)*diff;

if(initTorsion>10){
	cout<<"SugarVertex:: "<<initTorsion<<" .. "<<pos[4]<<" "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<" "<<endl;
	exit(-1);
}

	assert (fabs(initTorsion)<=Amplitude);
	assert (Amplitude>0.01);
}

double SugarVertex::getCurrentUP(){
	
	Vector3 p0 = ringAttached[0].front()->m_Position;
	Vector3 p1 = ringAttached[1].front()->m_Position;
	Vector3 p2 = ringAttached[2].front()->m_Position;
	Vector3 p3 = ringAttached[3].front()->m_Position;
	Vector3 p4 = ringAttached[4].front()->m_Position;
	Vector3 n = (p2-p1)/dists[1];
	Vector3 S1 = p1+n*( cos(CTK_PI-angles[2])*dists[2]+dists[1] );
	double S2off = n.dot(S1-p4);
	Vector3 S2 = p4+( n*( S2off ) );

	double p3side = cross( p3-S1,S2-S1 ).dot(n);
	return Signum(p3side);
}

/* Given the tau value fill the conf-array (assumed to have size 20) with all five coordinates 
and their reference frames. Returns true iff the chosen tau value (combined with the chosen amplitude) 
can result in closure of the sugar-ring. */
bool SugarVertex::getSugarConformation(double tau,Vector3* conf, double Am){//, bool usePosition2){
	
	double torsion = Am*cos(tau);
	double up = Signum(sin(tau));

	//Vector3 p0 = usePosition2?ringAttached[0].front()->Position2:ringAttached[0].front()->m_Position;
	//Vector3 p1 = usePosition2?ringAttached[1].front()->Position2:ringAttached[1].front()->m_Position;
	//Vector3 p4 = usePosition2?ringAttached[4].front()->Position2:ringAttached[4].front()->m_Position;
	Vector3 p0 = ringAttached[0].front()->m_Position;
	Vector3 p1 = ringAttached[1].front()->m_Position;
	Vector3 p4 = ringAttached[4].front()->m_Position;
	Vector3 p0x = p1-p0; normalize(p0x);
	Vector3 p0z = cross(p0x, p4-p0); normalize(p0z);
	Vector3 p0y = cross(p0z, p0x); 
	Vector3 p2 = p1+(cos(CTK_PI-angles[1])*dists[1])*p0x + (sin(CTK_PI-angles[1])*cos(torsion)*dists[1])*p0y + (sin(CTK_PI-angles[1])*sin(torsion)*dists[1])*p0z;

	Vector3 n = (p2-p1)/dists[1];
	//Vector3 S1 = p1+n*( cos(CTK_PI-angles[2])*dists[2]+dists[1] );
	Vector3 S1 = p2+n*(cos(CTK_PI-angles[2])*dists[2]);
	double S2off = n.dot(S1-p4);
	Vector3 S2 = p4 + n*S2off;
	double l = (S1-S2).length();
	double r1 = sin(CTK_PI-angles[2])*dists[2];
	double r2 = ( dists[3]*dists[3] - S2off*S2off );//Still need to take the square root
	if( r2<0 ) return false;
	r2 = sqrt(r2);
	//if(Am<0.1){
    //	log("debugRas")<<"Sphere["<<(p0).x<<","<<(p0).y<<","<<(p0).z<<",0.3,0.3,0.3,0.3]#C4"<<endl;
    //	log("debugRas")<<"Sphere["<<(p1).x<<","<<(p1).y<<","<<(p1).z<<",0.3,0.3,0.3,0.3]#C3"<<endl;
    //	log("debugRas")<<"Sphere["<<(p2).x<<","<<(p2).y<<","<<(p2).z<<",0.3,0.3,0.3,0.3]#C2"<<endl;
    //	log("debugRas")<<"Sphere["<<(p4).x<<","<<(p4).y<<","<<(p4).z<<",0.3,0.8,0.3,0.3]#O4"<<endl;
    //	log("debugRas")<<"Sphere["<<S1.x<<","<<S1.y<<","<<S1.z<<",0.2,0.3,0.3,0.8]#S1"<<endl;
    //	log("debugRas")<<"Sphere["<<S2.x<<","<<S2.y<<","<<S2.z<<",0.2,0.3,0.8,0.3]#S2"<<endl;
    //	log("debugRas")<<"Cylinder["<<S2.x<<","<<S2.y<<","<<S2.z<<","<<S2.x+0.05*n.x<<","<<S2.y+0.05*n.y<<","<<S2.z+0.05*n.z<<","<<r2<<",0.3,0.3,0.3]#S2"<<endl;
    //	log("debugRas")<<"Cylinder["<<S1.x<<","<<S1.y<<","<<S1.z<<","<<S1.x+0.05*n.x<<","<<S1.y+0.05*n.y<<","<<S1.z+0.05*n.z<<","<<r1<<",0.3,0.3,0.3]#S2"<<endl;
    //	log("debugRas")<<"r1 "<<r1<<endl;
    //	log("debugRas")<<"r2 "<<r2<<endl;
    //	log("debugRas")<<"l  "<<l<<endl;
    //	log("debugRas")<<"S2off "<<S2off<<endl;
    //	log("debugRas")<<"dists[4] "<<dists[4]<<endl;
    //	log("debugRas")<<"p0"<<p0<<endl;
    //	log("debugRas")<<"p1"<<p1<<endl;
    //	log("debugRas")<<"p2"<<p2<<endl;
    //	log("debugRas")<<"p4"<<p2<<endl;
	//}
	if(r1+r2<l) return false;
	if(conf==NULL) return true;

	double o = (r1+r2+l)/2;
	double A = sqrt(o*(o-r1)*(o-r2)*(o-l));
	double h = 2*A/l;
	double l1 = sqrt(r1*r1-h*h);
	Vector3 i = (S2-S1)/l;
	Vector3 j = cross(i,n);
	Vector3 p3 = S1+( i*l1 )+( j*up*h );


	Vector3 p3x = (p2-p3); normalize(p3x);
	Vector3 p3z = cross(p3x,p4-p3); normalize(p3z);
	Vector3 p3y = cross(p3z,p3x);
	Vector3 p2x = (p1-p2); normalize(p2x);
	Vector3 p2z = cross(p2x,p3-p2); normalize(p2z);
	Vector3 p2y = cross(p2z,p2x);
	Vector3 p1x = (p0-p1); normalize(p1x);
	Vector3 p1z = cross(p1x,p2-p1); normalize(p1z);
	Vector3 p1y = cross(p1z,p1x);
	Vector3 p4x = (p3-p4); normalize(p4x);
	Vector3 p4z = cross(p4x,p0-p4); normalize(p4z);
	Vector3 p4y = cross(p4z,p4x);

	conf[ 0] = p0;
	conf[ 1] = p0x;
	conf[ 2] = p0y;
	conf[ 3] = p0z;

	conf[ 4] = p1;
	conf[ 5] = p1x;
	conf[ 6] = p1y;
	conf[ 7] = p1z;

	conf[ 8] = p2;
	conf[ 9] = p2x;
	conf[10] = p2y;
	conf[11] = p2z;

	conf[12] = p3;
	conf[13] = p3x;
	conf[14] = p3y;
	conf[15] = p3z;

	conf[16] = p4;
	conf[17] = p4x;
	conf[18] = p4y;
	conf[19] = p4z;

	return true;
}

/*
	Updates the RigidbodyTransformation for edges going out from the sugar ring and for each adjacent Rigidbody. 
	It will also update the positions of atoms in ringAttached[2] and ringAttached[3]. 
*/
void SugarVertex::configurationToGlobalMatrix(RigidTransform* ms, double* m_f){//, bool usePosition2){
    //log("debugRas")<<"SugarVertex::configurationToGlobalMatrix( "<<relTau<<", ...) "<<endl;
	assert (fabs(initTorsion)<=Amplitude);
	double iniTau = acos(initTorsion/Amplitude)*initUp;
	double dTau = m_f[DOF_id];
	double curTau = iniTau+dTau;

	//Get the transformations for initial and new sugar-conformation
	Vector3 iniConf[20];
	Vector3 curConf[20];
    if(!getSugarConformation(iniTau,iniConf, Amplitude)){//, usePosition2)){
		cerr<<"SugarVertex::configurationToGlobalMatrix - Invalid sugar conformation .. adjust amplitude"<<endl; exit(-1); }
	if(!getSugarConformation(curTau,curConf, Amplitude)){//, usePosition2)){
		cerr<<"SugarVertex::configurationToGlobalMatrix - Invalid sugar conformation .. adjust amplitude"<<endl; exit(-1); }


	map<unsigned int,Edge*>::iterator edge_itr;
	RigidTransform localMat, m1, m2, m3, m4, m5,m6,m7;
	Vector3 vec;
	delayedUpdates.clear();

	//Update transformation for outgoing edges and positions of moving atoms within this rigidbody.
	for(int i=1;i<5;i++){
		m1.setIdentity();
		m2.setIdentity();
		m3.setIdentity();
		m4.setIdentity();
		m5.setIdentity();
		m6.setIdentity();
		m7.setIdentity();
		localMat.setIdentity();

		m1.setTranslate(curConf[i*4]);
		m2.R.set(curConf[i*4+1], curConf[i*4+2], curConf[i*4+3]);
		m3.R.set(iniConf[i*4+1], iniConf[i*4+2], iniConf[i*4+3]); m3.R.inplaceTranspose();//Orthonormal, so transpose==inverse
		m4.setTranslate(-1*iniConf[i*4]);
		localMat = m1*m2*m3*m4;

		for(list<Atom*>::iterator ait=ringAttached[i].begin(); ait!=ringAttached[i].end(); ait++){
			Atom* atm = *ait;
		//Atom* atm = *(ringAttached[i].begin());

			//If there is an outgoing edge from atm update its transformation
			//for (edge_itr=Edges.begin(); edge_itr != Edges.end(); ++edge_itr) {
			for (auto const& edge: edges) {
				//Edge* edge = (*edge_itr).second;
				//Edge* edge = *eit;
				if( edge->getBond()->Atom1==atm ){


					//if(!usePosition2) {
						m5.setTranslate(atm->m_Position);
						vec = edge->getBond()->Atom2->m_Position - edge->getBond()->Atom1->m_Position;
						m7.setTranslate(-1*atm->m_Position);
					//}else {
					//	m5.setTranslate(atm->Position2);
					//	vec = edge->getBond()->Atom2->Position2 - edge->getBond()->Atom1->Position2;
					//	m7.setTranslate(-1*atm->Position2);
					//}
					m6.setRotate(FindRotationMatrix(vec, -m_f[edge->DOF_id])); // !!! Since the FindRotationMatrix is for left hand, choose the negative of the angle

					ms[edge->DOF_id] = transformation * localMat * m5 * m6 * m7;
					edge->EndVertex->transformation = ms[edge->DOF_id];
				}
			}

			//Update the position of atm itself
			//if(!usePosition2)
				delayedUpdates[atm->getId()] = localMat*atm->m_Position;
			//else             	delayedUpdates[atm->getId()] = localMat*atm->Position2 ;
			//if(!usePosition2)	atm->m_Position.set( localMat*atm->m_Position );
			//else					atm->Position2.set( localMat*atm->Position2 );
		}
		//for(list<Atom*>::iterator ait=ringAttached[i].begin(); ait!=ringAttached[i].end(); ait++){
		//	Atom* atm = *(ait);
		//	if(!usePosition2) 	atm->m_Position.set(localMat*atm->m_Position);
		//	else					atm->Position2.set(localMat*atm->Position2);
		//}

	}
}

void SugarVertex::updateDelayed(){//bool usePosition2){
	for(vector<Atom*>::iterator ait=Rb_ptr->Atoms.begin(); ait!=Rb_ptr->Atoms.end(); ait++){
		Atom* atm = *ait;
		if(delayedUpdates.find(atm->getId())!=delayedUpdates.end()){
			//if(usePosition2) 	atm->Position2.set(delayedUpdates[atm->getId()]);
			//else
			atm->m_Position.set(delayedUpdates[atm->getId()]);
		}
	}
}

Vector3 SugarVertex::computeJacobianEntry(Edge* outEdge, double* m_f, Vector3& p){
	assert (fabs(initTorsion)<=Amplitude);
	double iniTau = acos(initTorsion/Amplitude)*initUp;
	double dTau = m_f[DOF_id];
	double curTau = iniTau+dTau;

	double delta = 0.0001;

	assert (!std::isnan(iniTau));

	
	Vector3 curConf[20];
	Vector3 newConf[20];
	if(!getSugarConformation(curTau,	  curConf, Amplitude)){ cerr<<"SugarVertex::computeJacobianEntry - Invalid sugar conformation .. adjust amplitude"<<endl; exit(-1); }
	if(!getSugarConformation(curTau+delta,newConf, Amplitude)){ cerr<<"SugarVertex::computeJacobianEntry - Invalid sugar conformation .. adjust amplitude"<<endl; exit(-1); }

	int confIdx = -1;
	for(int i=1;i<5;i++){
		if(std::find(ringAttached[i].begin(), ringAttached[i].end(), outEdge->getBond()->Atom1)!=ringAttached[i].end()){
			confIdx = i;break;
		}
	}
//	if(!(confIdx>1 && confIdx<5)) log("debugRas")<<confIdx<<endl;
//	assert(confIdx>1 && confIdx<5);

	//Vector3 a = p-curConf[confIdx*4];
	//Vector3 localCurP(a.dot(curConf[confIdx*4+1]), a.dot(curConf[confIdx*4+2]), a.dot(curConf[confIdx*4+3])); 
	//Vector3 localNewP = localCurP.x*newConf[confIdx*4+1]+localCurP.y*newConf[confIdx*4+2]+localCurP.z*newConf[confIdx*4+3] ;
	//Vector3 newP = localNewP+newConf[confIdx*4];
	//double dx = (newP.x-p.x)/delta;
	//double dy = (newP.y-p.y)/delta;
	//double dz = (newP.z-p.z)/delta;
	Vector3 ps = p-curConf[confIdx*4+0];
	Matrix3 r1(curConf[confIdx*4+1], curConf[confIdx*4+2], curConf[confIdx*4+3]); r1.inplaceTranspose(); 
	Matrix3 r2(newConf[confIdx*4+1], newConf[confIdx*4+2], newConf[confIdx*4+3]);
	Vector3 pe = r2*r1*ps + newConf[confIdx*4+0];
	Vector3 ret = (pe-p)/delta;
	//double* ret = new double[3];
	//ret[0] = dx;
	//ret[1] = dy;
	//ret[2] = dz;


    //log("debugRas")<<"Sphere["<<(curConf[ 0])[0]<<","<<(curConf[ 0])[1]<<","<<(curConf[ 0])[2]<<",0.3,0.3,0.3,0.3]#C1"<<endl;
    //log("debugRas")<<"Sphere["<<(curConf[ 4])[0]<<","<<(curConf[ 4])[1]<<","<<(curConf[ 4])[2]<<",0.3,0.3,0.3,0.3]#C2"<<endl;
    //log("debugRas")<<"Sphere["<<(curConf[ 8])[0]<<","<<(curConf[ 8])[1]<<","<<(curConf[ 8])[2]<<",0.3,0.3,0.3,0.3]#C3"<<endl;
    //log("debugRas")<<"Sphere["<<(curConf[12])[0]<<","<<(curConf[12])[1]<<","<<(curConf[12])[2]<<",0.3,0.3,0.3,0.3]#C4"<<endl;
    //log("debugRas")<<"Sphere["<<(curConf[16])[0]<<","<<(curConf[16])[1]<<","<<(curConf[16])[2]<<",0.3,0.8,0.3,0.3]#O4"<<endl;

    //log("debugRas")<<"Sphere["<<(newConf[ 0])[0]<<","<<(newConf[ 0])[1]<<","<<(newConf[ 0])[2]<<",0.3,0.3,0.3,0.3]#C1"<<endl;
    //log("debugRas")<<"Sphere["<<(newConf[ 4])[0]<<","<<(newConf[ 4])[1]<<","<<(newConf[ 4])[2]<<",0.3,0.3,0.3,0.3]#C2"<<endl;
    //log("debugRas")<<"Sphere["<<(newConf[ 8])[0]<<","<<(newConf[ 8])[1]<<","<<(newConf[ 8])[2]<<",0.3,0.3,0.3,0.3]#C3"<<endl;
    //log("debugRas")<<"Sphere["<<(newConf[12])[0]<<","<<(newConf[12])[1]<<","<<(newConf[12])[2]<<",0.3,0.3,0.3,0.3]#C4"<<endl;
    //log("debugRas")<<"Sphere["<<(newConf[16])[0]<<","<<(newConf[16])[1]<<","<<(newConf[16])[2]<<",0.3,0.8,0.3,0.3]#O4"<<endl;

    //log("debugRas")<<"Sphere["<<p[0]<<","<<p[1]<<","<<p[2]<<",0.3,0.8,0.3,0.3]#O4"<<endl;
    //log("debugRas")<<"Sphere["<<newP[0]<<","<<newP[1]<<","<<newP[2]<<",0.3,0.3,0.8,0.3]#O4"<<endl;


	return ret;
}

