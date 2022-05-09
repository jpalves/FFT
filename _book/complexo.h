//Jo�o Pedro Alves 2004(c)
//correio para bugs jpalves@mail.isec.pt

#ifndef __cplusplus
	#error Tem de compilar em C++ para usar o tipo complexo.
#endif

#define __COMPLEXO_H

//erros inerentes � truncatura de dizimas infinitas
#define ERRO  1e-15

#if !defined(__MATH_H)
	#include <cmath>
#endif

#if !defined(__STDLIB_H)
	#include <stdlib.h>
#endif

#if !defined(__IOSTREAM_H)
	#include <iostream>
#endif

#include <string>

#define negativo(in) (in<0?1:0)


/*definir pi*/
#ifndef M_PI 
        #define M_PI		3.14159265358979323846
#endif


//dispon�vel em string.h com o nome strlen
int tamanhoStr(char *str){
	char *ptr_str;

	for(ptr_str=str;*ptr_str!='\0';ptr_str++);
	return (int)(ptr_str - str);
}
//existe em ctype.h com o mesmo nome
char tolower(char in){
	return in|(char)0x20;
}
//insere um caractere numa string na posi��o desejada 
void insere(char in,char *str,int posicao){
	char *ptr_str=str+tamanhoStr(str);

	*(ptr_str+1)='\0';
	for(;(str + posicao) - ptr_str;ptr_str--) *ptr_str=*(ptr_str -1);
	*ptr_str=in;
}

//---------------------------------------------------------------------------
class complexo{
	private:
		//campos
		double _real,imaginario;

		//fun��es friend n�o membros
		friend complexo cos(complexo z);
		friend complexo sin(complexo z);
		friend complexo exp(complexo z);

		//operadores friend
		friend complexo operator /(double real,complexo z);
		friend complexo operator /(complexo z,double real);
		friend complexo operator *(double real,complexo z);
		friend complexo operator *(complexo z,double real);
		friend complexo operator +(double real,complexo z);
		friend complexo operator -(double real,complexo z);
		friend std::ostream &operator<<(std::ostream &stream,complexo z);
		friend std::istream &operator>>(std::istream &stream,complexo &z);

	public:
		//construtores
		complexo(double _real,double imaginario=0){
			this->_real=_real;
			this->imaginario=imaginario;
		}
		complexo(){
			_real=0;imaginario=0;
		}
		complexo(const complexo &z){
			this->_real=z._real;
			this->imaginario = z.imaginario;
		}
		//m�todos
		inline double   modulo(){return sqrtl(powl(_real,2)+powl(imaginario,2));}
		inline double   abs()   {return sqrtl(powl(_real,2)+powl(imaginario,2));}
		inline const double abs() const {return sqrtl(powl(_real,2)+powl(imaginario,2));}
		       double   arg();
		inline double   argGrau() {return 180*arg()/M_PI;}
		inline double   &real()   {return _real;}
		inline double   &imag()   {return imaginario;}
		inline complexo      conj()    {return complexo(_real,-imaginario);}
		int                  quadrante();
		void                 setPolar(double modulo,double angulo){
			if(fabsl(cosl(angulo))>ERRO)_real=modulo*cosl(angulo);
			else _real=0;
			if(fabsl(sinl(angulo))>ERRO) imaginario=modulo*sinl(angulo);
			else imaginario=0;
		}
		void                 setPolarGrau(double modulo,double angulo){
			if(fabsl(cosl(M_PI*angulo/180))>ERRO) _real=modulo*cosl(M_PI*angulo/180);
			else _real=0;
			if(fabsl(sinl(M_PI*angulo/180))>ERRO) imaginario=modulo*sinl(M_PI*angulo/180);
			else imaginario=0;
		}
		void                 setRectangular(double re,double im){
			_real=re;
			imaginario=im;
		}
		//operadores reescritos
		complexo operator  /(complexo z);
		complexo operator  *(complexo z);
		complexo operator  +(complexo z);
		complexo operator  +(); //un�rio
		complexo operator  +(double real);
		complexo operator  -(complexo z);
		complexo operator  -(); //un�rio
		complexo operator  -(double real);
		bool     operator ==(complexo z);
		bool     operator !=(complexo z);
		complexo operator  =(complexo z);
		complexo operator +=(complexo z);
		complexo operator -=(complexo z);
		complexo operator *=(complexo z);
		complexo operator  =(double real);
		complexo operator /=(complexo z);
		complexo operator  ^(double real);
		bool     operator  !(); //not
		inline   operator void *(){return _real||imaginario?this:NULL;}
		inline   bool operator <(complexo z) const { return (this->abs() < z.abs());}
};
//---------------------------------------------------------------------------
inline complexo complexo::operator +(){
	return *this;
}

inline complexo complexo::operator +(complexo z){
	return complexo(_real+z._real,imaginario+z.imaginario);
}

inline complexo complexo::operator +(double real){
	return complexo(_real+real,imaginario);
}

inline complexo complexo::operator -(complexo z){
	return complexo(_real-z._real,imaginario-z.imaginario);
}

inline complexo complexo::operator -(){
	return complexo(-_real,-imaginario);
}

inline complexo complexo::operator -(double real){
	return complexo(_real-real,imaginario);
}
//divis�o de complexos
inline complexo complexo::operator /(complexo z){
	double modulo2=(z._real*z._real+z.imaginario*z.imaginario);

	return complexo((_real*z._real+imaginario*z.imaginario)/modulo2,(imaginario*z._real-_real*z.imaginario)/modulo2);
}
//multiplica��o de complexos
inline complexo complexo::operator *(complexo z){
	return complexo(_real*z._real-imaginario*z.imaginario,_real*z.imaginario+z._real*imaginario);
}
//atribui��o
complexo complexo::operator =(complexo z){

	_real=z._real;
	imaginario=z.imaginario;
	return *this;
}

complexo complexo::operator =(double real){

	_real=real;
	imaginario=0;
	return *this;
}

complexo complexo::operator +=(complexo z){

	 _real+=z._real;
	 imaginario+=z.imaginario;
	 return *this;
}

complexo complexo::operator -=(complexo z){

	 _real-=z._real;
	 imaginario-=z.imaginario;
	 return *this;
}

complexo complexo::operator *=(complexo z){

	*this=*this*z;
	return *this;
}

complexo complexo::operator /=(complexo z){

	*this=*this/z;
	return *this;
}

inline bool complexo::operator ==(complexo z){
	return !(_real-z._real)&&!(imaginario-z.imaginario);
}

inline bool complexo::operator !=(complexo z){
	return _real-z._real||imaginario-z.imaginario;
}

inline bool complexo::operator !(){
	return !(_real||imaginario);
}

complexo complexo::operator  ^(double expoente){
	complexo temp;

	temp.setPolar(powl(this->modulo(),expoente),this->arg()*expoente);
	return temp;
}
//---------------------------------------------------------------------------
//operador para ostream (ficheiro ou ecran)

std::ostream &operator <<(std::ostream &stream,complexo z){
        		
	if(!(z.imag()-1)&&!z.real()){ stream <<'i'; return stream;}
	if(!(z.imag()+1)&&!z.real()){ stream <<"-i"; return stream;}
	if(!(z.imag()-1))           { stream <<z.real()<<"+i"; return stream;}
	if(!(z.imag()+1))           { stream <<z.real()<<"-i"; return stream;}
	if(!z.imag())               { stream << z.real(); return stream;}
	if(!z.real())               { stream <<z.imag()<<"i"; return stream;}
	if(z.imag() < 0)stream << z.real()<<z.imag()<<'i';
	else stream << z.real()<<'+'<<z.imag()<<'i';
	return stream;
}
//operador para istream (ficheiro ou teclado)
std::istream &operator >>(std::istream &stream,complexo &z){
	char str[256]="",*ptr=str,imag[128]="";

	stream >> str;
	z.imaginario=0;

	for(int i=0;str[i];i++)

	if(tolower(str[i])=='i'||tolower(str[i])=='j'){
		*(ptr+i)=(char)0x20;
		for(;&ptr[i] != str;ptr--){
			insere(*(ptr+i-1),imag,0);
			ptr[i-1]=(char)0x20;
			if(( *imag == '+')||( *imag == '-')||(&ptr[i-1] == str)) break;
		}
		ptr=str;
		if(tamanhoStr(imag)>1){z.imaginario=(double)strtod(imag, NULL);} //atof em android dá bronca
		else{
			if(*imag!='\0'&&(*imag != '+')&&(*imag != '-')) z.imaginario=strtod(imag, NULL);
			if((*imag=='-')) z.imaginario=-1;
			if((*imag=='+')) z.imaginario=1;
			if(*imag =='\0') z.imaginario=1;
		}
		*imag='\0';
	}
	z._real=(double)strtod(str, NULL);
	return stream;
}
//---------------------------------------------------------------------------
//m�todos
int complexo::quadrante(){

	if(!_real&&!imaginario)           return 0; //refer�ncia
	if(!_real&&!negativo(imaginario)) return 12;//fronteira entre o 1� e 2� q
	if(!_real&&negativo(imaginario))  return 34;
	if(!imaginario&&negativo(_real))  return 23;
	if(!imaginario&&!negativo(_real)) return 14;
	if(!negativo(_real)&&!negativo(imaginario)) return 1;
	if(negativo(_real)&&!negativo(imaginario))  return 2;
	if(negativo(_real)&&negativo(imaginario))   return 3;
	if(!negativo(_real)&&negativo(imaginario))  return 4;
	return 0;
}

double complexo::arg(){
	double aux;

	switch(quadrante()){
		case  0: aux= 0;                            break;
		case 12: aux= M_PI/2;                       break;
		case 23: aux= M_PI;                         break;
		case 34: aux=-M_PI/2;                       break;
		case  2: aux= M_PI+atanl(imaginario/_real); break;
		case  3: aux= atanl(imaginario/_real)-M_PI; break;
		default: aux= atanl(imaginario/_real);
	}
	return aux;
}

//---------------------------------------------------------------------------
//fun��es friend
//cosseno de um complexo
/*
	  jz   -jz        j(a+bj)       -j(a+bj)      -b                     b
	(e  + e    )/2= (e      +      e        )/2= (e *(cos a + j sen a)+ e (cos a + j sen a))/2
									    |
									    equivalente a cos -a
	basta agora separar a parte real e a parte imagin�ria
*/
inline complexo cos(complexo z){
	return complexo(expl(-z.imaginario)*cosl(z._real)/2.0+expl(z.imaginario)*cosl(-z._real)/2.0,expl(-z.imaginario)*sinl(z._real)/2.0+expl(z.imaginario)*sinl(-z._real)/2.0);
}
//seno de um complexo
inline complexo sin(complexo z){
	return complexo(expl(-z.imaginario)*sinl(z._real)/2.0-expl(z.imaginario)*sinl(-z._real)/2.0,-expl(-z.imaginario)*cosl(z._real)/2.0+expl(z.imaginario)*cosl(-z._real)/2.0);
}
//expon�ncial complexa
// ver formula de Euler
inline complexo exp(complexo z){
	return complexo(expl(z._real)*cosl(z.imaginario),expl(z._real)*sinl(z.imaginario));
}
//---------------------------------------------------------------------------
//divis�o de um  n� real por um complexo
complexo operator /(double real,complexo z){
	double z_conjz=z._real*z._real+z.imaginario*z.imaginario;//powl(z.modulo(),2);

	return complexo(real*z._real/z_conjz,-real*z.imaginario/z_conjz);
}

inline complexo operator /(complexo z,double real){
	return complexo(z._real/real,z.imaginario/real);
}

inline complexo operator *(complexo z,double real){
	return complexo(z._real*real,z.imaginario*real);
}

inline complexo operator *(double real,complexo z){
	return complexo(z._real*real,z.imaginario*real);
}

inline complexo operator +(double real,complexo z){
	return complexo(real+z._real,z.imaginario);
}

inline complexo operator -(double real,complexo z){
	return complexo(real - z._real,- z.imaginario);
}
//---------------------------------------------------------------------------
//fun��es normais overloded (reescritas)
//ramo principal do logar�tmo
complexo log(complexo z){
	 complexo temp;

	 temp.setRectangular(log(z.modulo()),z.arg());//z.arg()+2k*pi k pretence n naturais
	 return temp;
}

complexo tan(complexo z){
	complexo temp=sin(z)/cos(z);

	return temp;
}

complexo sh(complexo z){
	complexo j(0,1);// 1 cis(pi/2) = i

	return sin((-j)*z)*j;
}

complexo csh(complexo z){      //cosseno hiperb�lico
	complexo j(0,1);

	return cos((-j)*z);
}

complexo pow(complexo z,long double expoente){
	 complexo temp;

	 temp.setPolar(powl(z.modulo(),expoente),z.arg()*expoente);
	 return temp;
}

//primeira raiz
//a segunda raiz pode ser obtida somando ao argumento +2K*pi/n  em que n=2 e K=1.
complexo sqrt(complexo z){
	 return pow(z,(double)1/2);
}

complexo atan(complexo z){
	complexo j(0,1);

	return 0.5*(-j)*log((1 +j*z)/(1 + (-j)*z));
}

complexo ash(complexo z){
	return log(z+ sqrt(pow(z,2) + 1));
}

complexo acsh(complexo z){
	return log(z +sqrt(pow(z,2) - 1));
}

complexo asin(complexo z){
	complexo j(0,1);

	return (-j)*log(j*z +sqrt(-z*z + 1));
}
//c�digo exprimental (argumento do cosseno hiperb�lico)
complexo acos(complexo z){
       complexo j(0,1);

       if(!(z.quadrante()-2)||!(z.quadrante()-4)) return j*log(z + sqrt(z*z - 1));
       return (-j)*log(z + sqrt(z*z - 1));
}
