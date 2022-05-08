#include <iostream>
#include <cassert>
#include <limits>

//
template <typename T> class matriz{
	int m,n;
	T **v;
	//inline void swap(int &a,int &b){int aux=a;a=b;b=aux;}
	template <typename U> friend matriz<U> operator /(const matriz<U> &in,const U &z);
        template <typename U> friend matriz<U> operator *(const U &z, const matriz<U> &in);	
public:
	matriz():m(0),n(0),v(NULL){}
	matriz(int m,int n);
	matriz(int m,int n,const T *a);
	matriz(int m,int n,const T &a);
	matriz(const matriz &rhs);	
	~matriz();
	matriz<T>  operator*(const matriz<T> &in);
	matriz<T>  operator+(const matriz<T> &in);
    	matriz<T>  operator-(const matriz<T> &in);
	matriz<T>  operator*(const T s);
	inline T*  operator[](const int i){return v[i];}
	matriz<T> &operator=(const matriz<T> &in);
	T  *asVector()      {return v[0];}
	int linhas()  const {return m;}
	int colunas() const {return n;} 
	int numel()   const {return m*n;}
    	T det();
	matriz<T> inv();
    	matriz<T> coff(unsigned i, unsigned j) const;
    	matriz<T> Inv();
    	matriz<T> transposta();
    	void trocaLinhas(int i,int j);
    	void trocaColunas(int i,int j);
    	matriz<T> eye();	
};

template <typename T> matriz<T>::matriz(int m,int n):m(m),n(n),v(n > 0 ?new T *[m]: NULL){
	int nEl = m*n;
		
	if(v)
		v[0] = nEl > 0 ? new T[nEl] : NULL;
	for(int i=1;i < m ;i++) v[i] = v[i-1] + n;
	//std::cout <<m<<','<<n<<std::endl;
}

template <typename T> matriz<T>::matriz(int m,int n,const T *a):m(m),n(n),v(m > 0 ?new T *[m]: NULL){
	int nEl = m*n;
		
	if(v)
		v[0] = nEl > 0 ? new T[nEl] : NULL;
	for(int i=1;i < m ;i++)  v[i] = v[i-1] + n;
	for(int i=0;i < nEl;i++) v[0][i] = *a++;
}

template <typename T> matriz<T>::matriz(int m,int n,const T &a):m(m),n(n),v(m > 0 ?new T *[m]: NULL){
	int nEl = m*n;
		
	if(v)
		v[0] = nEl > 0 ? new T[nEl] : NULL;
	for(int i=1;i < m ;i++)  v[i] = v[i-1] + n;
	for(int i=0;i < nEl;i++) v[0][i] = a;
}

template <typename T> matriz<T>::matriz(const matriz &rhs):m(rhs.m),n(rhs.n),v(m > 0 ?new T *[m]: NULL){
	int nEl = m*n;
	
	if(v)
		v[0] = nEl > 0 ? new T[nEl] : NULL;
	for(int i=1;i < m ;i++)  v[i] = v[i-1] + n;
	for(int i=0;i < nEl;i++) v[0][i] = rhs.v[0][i];

}

template <typename T> matriz<T>::~matriz(){
	if (v != NULL) {
		delete[] (v[0]);
		delete[] (v);
	}
}

template <typename T> matriz<T> matriz<T>::operator *(const matriz<T> &in){
	if(!(n-in.m)){	
		matriz<T> temp(m,in.n,(T)0);
		for(int i=0;i < m;i++)
			for(int j=0;j < in.n;j++)
				for(int k=0;k < n;k++){
					temp.v[i][j] += v[i][k]*in.v[k][j];
				}
		return temp;
	} else return matriz<T>();
	
	
}

template <typename T> matriz<T> matriz<T>::operator +(const matriz<T> &in){
	matriz<T> temp(m,n);
	
	//falta proteger
	for(int i = 0; i < m*n;i++)
		temp.v[0][i] = v[0][i] + in.v[0][i];
	return temp;
}

template <typename T> matriz<T> matriz<T>::operator -(const matriz<T> &in){
	matriz<T> temp(m,n);
	
	//falta proteger
	for(int i = 0; i < m*n;i++)
		temp.v[0][i] = v[0][i] - in.v[0][i];
	return temp;
}

//código manhoso
template <typename T>matriz<T> &matriz<T>::operator =(const matriz &in){
	//std::cout <<"entrei"<<std::endl;
	if(this != &in){
		int nEl;
		//std::cout <<"entrei"<<std::endl;
		//std::cout <<m<<n<<in.m<<in.n<<std::endl;
		if (n != in.n || m != in.m) {
			//std::cout <<"entrei"<<std::endl;
			if (v != NULL) {
				delete[] (v[0]);
				delete[] (v);
			}
			n=in.n;
			m=in.m;
			v = m > 0 ? new T*[m] : NULL;
			nEl = m*n;
			if(v) //o código de alocção e desalocação vai para uma função privada
				v[0] = nEl > 0 ? new T[nEl] : NULL;
			for(int i=1;i < m ;i++) v[i] = v[i-1] + n;
		}
		for(int i=0;i < m*n;i++) v[0][i] = in.v[0][i]; //verificar este código
		//std::cout <<"sai"<<std::endl;
	}
	return *this;
}

template <typename U> matriz<U> operator /(const matriz<U> &in,const U &z){
	matriz<U> temp = in;
	
	for(int i = 0;i < temp.m*temp.n;i++) temp.v[0][i] /=z;
	return temp;
}


template <typename T> matriz<T> matriz<T>::operator*(const T s){
	matriz<T> temp = *this;
	
	for(int i=0; i < m; i++)
		for(int j=0; j < n; j++)
			temp[i][j] *= s; 
        return temp;
}


template <typename U> matriz<U> operator *(const U &z,const matriz<U> &in){
	matriz<U> temp = in;
	
	for(int i = 0;i < temp.m*temp.n;i++) temp.v[0][i] *=z;
	return temp;
}

template <typename T> T matriz<T>::det(){
        T aux = 1,tempC;
	int i,j,k,l;

	if(m-n){std::cout <<"A matriz não é quadrada"<<std::endl; return 0;}
	else{//Algoritmo de Gauss para condensação de matrizes
		matriz<T> temp(m,n);
		temp=*this;

		for(k=0;k<temp.m-1;k++){//pivot
			if(!temp[k][k]){ //troca de linhas caso o pivot seja nulo
				for(l=k+1;l<temp.m;l++)
					if(temp[l][k] != (T) 0)break;

				for(j=k;j<temp.n;j++){
					tempC=temp[k][j];
					temp[k][j]=-temp[l][j];
					temp[l][j]=tempC;
				}
			}
			
			for(i=k+1;i<temp.m;i++)
				if(temp[k][k]){
					for(j=k+1;j<temp.n;j++)
						temp[i][j] -= (temp[i][k]*temp[k][j])/temp[k][k];
				}
		}
	       //determinante é o produto dos elementos da diagonal principal
	       for(i=0;i<temp.m;i++)aux *= temp[i][i];
	}
	return aux;
}

template <typename T> matriz<T> matriz<T>::inv(){
	matriz<T> temp(m,n),iter(m,n);
        matriz <bool> premut(m,n);
	int i,j,k,l;
        

	if(m - m){std::cout <<"A matriz não é quadrada"<<std::endl;}
	else{
		//if(!det()){std::cout <<"A matriz é singular"<<std::endl; return temp;}
		assert(det() != (T)0);
                temp=*this;
		for(k=0;k < m; k++){
			if(!v[k][k]){
				for(l=k+1;;l++){
					if(!(l-m))l=0;
					if((l-k)&&v[l][k])break;
				}
			     temp.trocaLinhas(k,l);
			     premut[k][k]=true;
			     premut[k][l]=true;
			}

			iter[k][k]= -1/temp[k][k];
			for(i=0;i < m;i++)
				for(j=0;j < n;j++){
                             		if( (i-k)&& (j-k)) iter[i][j]  = temp[i][j] + temp[k][j]*temp[i][k]*iter[k][k];
                               	if( (i-k)&&!(j-k)) iter[i][k]  = temp[i][k]*iter[k][k];		
					if(!(i-k)&& (j-k)) iter[k][j]  = temp[k][j]*iter[k][k];                          
				}
			for(i=0;i < m;i++) for(j=0;j < n;j++){
							temp[i][j]=iter[i][j];
						   }
		}
		for(i=0;i<temp.m;i++)
			for(int j=0;j<temp.n;j++)
				temp[i][j]=-temp[i][j];


		for(i=temp.m-1;i+1;i--)
			 if(premut[i][i]){
				for(k=0;;k++){
					if(!(k-m))k=0;
					if((k-i)&&premut[i][k])break;
				}
				temp.trocaColunas(i,k);
			 }
      }
      return temp;
}

template <typename T> void matriz<T>::trocaLinhas(int i, int j){
	T aux;

	for(int k=0;k < n;k++){
		aux=this->v[i][k];
		this->v[i][k]=this->v[j][k];
		this->v[j][k]=aux;
	}
}

template <typename T> void matriz<T>::trocaColunas(int i, int j){
	T aux;

	for(int k=0;k < m;k++){
		aux=this->v[k][i];
		this->v[k][i]=this->v[k][j];
		this->v[k][j]=aux;
	}
}

//
template <typename T> matriz<T> matriz<T>::coff(unsigned i, unsigned j) const {
	if (n==0)
		throw std::logic_error("Matriz vazia");

	matriz<T> y(n-1,m-1);
	
	unsigned k_c=0;
	for (unsigned k_x=0; k_x< n;k_x++){
		if (k_x==i)
			continue;
		unsigned j_c=0;
		for (unsigned j_x=0 ; j_x<n; j_x++){
			if (j_x==j)
				continue;
			y.v[k_c][j_c]=v[k_x][j_x];
			j_c++;
		}
		k_c++;
	}
	return y;	
}

//
template <typename T> matriz<T> matriz<T>::Inv(){
	double det_x=det();
	if(abs(det_x)<std::numeric_limits<double>::epsilon())
		throw std::logic_error("Matriz Singular impossível inverter");
	
	matriz<T> y(n,m);
	
	signed int d=1;

	for (unsigned i=0; i<n; i++)
		for (unsigned j=0; j<n; j++){
			y.v[j][i]=d*coff(i,j).det()/det_x;
			d=-d;
		}

	return y;
	
}

template <typename T> matriz<T> matriz<T>::transposta(){
	matriz<T> y(m,n);

        for (unsigned i=0; i< n; i++)
		for (unsigned j=0; j< m; j++)
			y.v[j][i]=v[i][j];
	return y;
}

template <typename T> matriz<T> matriz<T>::eye(){
	matriz<T> y(m,n);

        for (unsigned i=0; i< n; i++)
		y.v[i][i]=1;
	return y;
}  

