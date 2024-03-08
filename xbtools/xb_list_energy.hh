///This class is based on a suggestion by Tudi. Thanks, Tudi :)
#ifndef XB_LIST_ENERGY
#define XB_LIST_ENERGY

#define GAMMA_MAX_THRESHOLD 25.0
#define GAMMA_MIN_THRESHOLD 0.6
#define PROTON_THRESHOLD 25.0
#define T_MIN -80.0
#define T_MAX -45.0

#include <iostream>
#include <list>
#include <string>
#include <math.h>

//for making lists of crystal (hit)s, sortable by energy
class xble
{

	public:
		xble() { reset(); }

		//xble(int n, float m, float o, float q)
		xble(int n, float m, float o, float q, int mul)
		{
			this->set(n,m,o,q,mul);
			//this->set(n,m,o,q);
		}

		xble(const xble &r)
		{
			//this->set(r._number, r._energy_p, r._energy_g, r._time);
			this->set(r._number, r._energy_p, r._energy_g, r._time, r._mul);
		}

	public:
		~xble() { }

	public:
		void reset()
		{
			this->_number = -1;
			this->_energy_p = -1;//(double)nan;
			this->_energy_g = -1;//(double)nan;
			this->_time = -1;//(double)nan;
			this->_mul = 0;
		}

	public:
		//void set(int n, float m, float o, float q)
		void set(int n, float m, float o, float q, int mul)
		{
			this->_number = n;
			this->_energy_p = m;
			this->_energy_g = o;
			this->_time = q;
			this->_mul = mul;
		}

	public:
		bool operator==(const xble &r) const
		{
			if(this->_number != r._number) return false;
			if(this->_energy_p != r._energy_p) return false;
			if(this->_energy_g != r._energy_g) return false;
			if(this->_time != r._time) return false;
			return true;
		}

		bool operator < (const xble &r) const
		{
			//if(this->_energy_p < r._energy_p) return 1;

			//This is a manually 'crafted' sorting definition, matching the needs for sorting energies in the two XB readouts:

			//First, sort by gamma energy entry. If that's inf, then by proton energy entry. Later need to throw out events where gamma was inf, but proton nan, maybe also where there was proton, but nan in gamma...

			//first check: sort by gamma energy
			if(this->_energy_g > r._energy_g) return true;
			
			if (isinf(this->_energy_g) && isinf(r._energy_g))
			{
				//if gamma in overflow, check protons
				if (this->_energy_p > r._energy_p) return true;
			}
			
			//if one gamma entry is nan, it shall be treated as smaller
			if(!isnan(this->_energy_g)&&isnan(r._energy_g)) return true;
			
			if(isnan(this->_energy_g)&&isnan(r._energy_g))
			{
				//if both gamma entries are nan, go for proton entry
				if (this->_energy_p > r._energy_p) return true;
				//if(!isnan(this->_energy_p)&&isnan(r._energy_p)) return 1;
			}
			return false;
		}


		xble &operator=(const xble &r)
		{
			//this->set(r._number, r._energy_p, r._energy_g, r._time);
			this->set(r._number, r._energy_p, r._energy_g, r._time, r._mul);

			return *this;
		}

	public:
		void print(std::ostream &oss=std::cout) const
			//void print(std::ostream &oss=std::cout)
		{
			oss << "Crystal " << this->_number << ":\t" << this->_energy_p << "\t"<< this->_energy_g << "\t" << this->_time << "\t" << this->_mul <<std::endl;
			//oss << "Crystal " << this->_number << ":\t" << this->_energy_p << ":\t" << this->_energy_g << ":\t" << this->_mul << std::endl;
		}

	public:
		//function for add-back routine
		//bool is_proton_4pi() const
		//{
		//	if ((_energy_p>PROTON_THRESHOLD && (_energy_g>GAMMA_MAX_THRESHOLD || isinf(_energy_g) || isnan(_energy_g))) || 
		//			(_energy_p<PROTON_THRESHOLD && isinf(_energy_g)) ||
		//		(isnan(_energy_p) && isinf(_energy_g))) //included overflow in back hemisphere!!!!!!!
		//		return true;
		//	else
		//		return false;
		//}

		bool is_proton_4pi() const
		{
			if((_energy_p>PROTON_THRESHOLD && (_energy_g>GAMMA_MAX_THRESHOLD || isinf(_energy_g))) || 
					(_energy_p<PROTON_THRESHOLD && isinf(_energy_g)) ||
				(isnan(_energy_p) && isinf(_energy_g))) //included overflow in back hemisphere!!!!!!!
				return true;
			else
				return false;
		}


	public:
		//function for add-back routine
		bool is_proton_2pi() const
		{
			if (_energy_p>PROTON_THRESHOLD && !isinf(_energy_p) && (_energy_g>GAMMA_MAX_THRESHOLD || isinf(_energy_g) || isnan(_energy_g)))
				return true;
			else
				return false;
		}

		//function for add-back routine
		bool is_proton() const
		{
			if (_energy_p>0.0 && (_energy_g>GAMMA_MAX_THRESHOLD || isinf(_energy_g) || isnan(_energy_g)))
				return true;
			else
				return false;
		}




	public:
		//function for add back routine
		bool is_gamma() const
		{
			if (_energy_g<GAMMA_MAX_THRESHOLD || (_energy_g>GAMMA_MAX_THRESHOLD && !isinf(_energy_g) &&  isnan(_energy_p)))
				return true;
			else
				return false;
		}
	
	public:
		//function for add back routine
		bool is_gamma(Float_t t0) const
		{
			if ((_energy_g<GAMMA_MAX_THRESHOLD || (_energy_g>GAMMA_MAX_THRESHOLD && !isinf(_energy_g) &&  isnan(_energy_p))) 
				&& (_time-t0)>T_MIN && (_time-t0)<T_MAX) //time window for "good" gammas
				return true;
			else
				return false;
		}



	public:
		//private:
		int	_number;
		int	_mul;
		float 	_energy_p;
		float 	_energy_g;
		float 	_time;

};

//-------------------------------------------

void print_all(const std::list<xble> &rl);

#endif //XB_LIST_ENERGY
