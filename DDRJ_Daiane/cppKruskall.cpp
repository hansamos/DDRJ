#include <Rcpp.h>
using namespace Rcpp;

//Informações sobre a função kruskalcpp
//kruskalcpp é uma função codificada em C++ e tem como objetivo realizar o teste Kruskal Wallis presente na linguagem R.
//A entrada da função kruskalcpp é uma matriz.
//Sua saída é um vetor com dimensão equivalente a diferença entre o número de colunas da matriz de entrada menos 1 (nº de colunas - 1).

// [[Rcpp::export]]
NumericVector kruskalcpp( NumericMatrix dados ){
  
  //"n_linha" equivale ao número de linhas presentes na matrix de entrada.
  //"n_linha" não é do tipo inteiro, pois será usada para calculos envolvendo variáveis do tipo float.
  float n_linha;
  n_linha = dados.nrow();
  
  //"n_col" corresponde ao número de colunas presentes na matrix de entrada.
  int n_col;
  n_col = dados.ncol();
  
  //Bloco encarregado por computar o rank dos elementos que compõem a primeira coluna da matriz de entrada.
  //Devolverá um vetor nomeado como "rank" composto pela colocação crescente de cada elemento da primeira coluna do vetor de entrada.
  //Exemplo:
  //1ª coluna do vetor_entrada : -4.4 5.6 -28.9 -8.7
  //rank : 3 4 1 2
  NumericVector rank( n_linha );
  for ( int i = 0; i < n_linha; i++){
    
    int rank_temporario = 1;
    
    for ( int j = 0; j < i; j++){
      
      if ( dados( i , 0 ) > dados( j, 0)){
        rank_temporario++;
      }
      else{
        rank[ j ]++;
      }
      
    }
    rank[ i ] = rank_temporario;
  }
  
  //Vetor "kw" armazena os valores do teste de Kruskal Wallis, sua dimensão equivale a diferença do número de colunas do vetor de entrada menos 1.
  //"tam_kw" equivale ao tamanho do vetor kw.
  int tam_kw = n_col - 1;
  NumericVector kw( tam_kw );
  
  //Vetor "freq_ale" contém a frequência de cada tipo de alelo presente na matrix de entrada.
  //Possui tamanho igual a 3, pois os alelos assumem valores iguais a -1, 0 e 1.
  //Onde o primeiro elemento de "freq_ale" corresponde a frequência absoluta do alelo "-1", o segundo a "0" e o terceiro "1" na matrix de entrada.
  NumericVector freq_ale(3);
  
  //O vetor "soma_rank" contém a soma total do rank para cada tipo de alelo.
  //Possui tamanho igual a 3.
  NumericVector soma_rank(3);
  
  //"razao_rank_freqale" armazena a razao entre o quadrado de soma_rank e freq_ale.
  NumericVector razao_rank_freale(3);
  
  //"total" é a soma das razões entre o quadrado de "soma_rank" e "freq_ale" para cada tipo de alelo.
  float total( n_col - 1);
  
  //Laço responsável por inserir os valores que formarão os elementos do vetor kw, tal estrura de dado armazena valores do teste de Kruskal Wallis.
  for ( int j = 1; j < n_col; j++){
    
    //A cada vez que o programa varrer uma nova coluna,inicialmente, zero será atribuído a variável total
    total = 0;
    
    //Laço iguala a zero todos os elementos dos vetores "freq_ale" e "soma_rank", a cada vez que o programa varrer uma nova coluna da matriz de entrada.
    for ( int m = 0; m < 3; m++ ){
      
      freq_ale( m ) = 0;
      
      soma_rank( m ) = 0;
      
    }
    
    //Laço percorre todas as linhas da matrix de entrada.
    for ( int i= 0; i < n_linha; i++){
      
      //Laço checa qual é o tipo de alelo, computa sua frequência e soma o rank referente a cada tipo de alelo.
      for ( int k = -1; k < 2; k++ ){
        
        if ( dados( i , j ) == k ){
          freq_ale( k + 1 ) += 1;
          soma_rank( k + 1 ) += rank( i ); 
        }
        
      }
      
      for ( int n = 0; n < 3; n++ ){ 
        
        if ( soma_rank( n ) != 0){
          razao_rank_freale( n ) =  ( soma_rank( n ) * soma_rank( n ) ) / freq_ale( n );
        }
        else{
          razao_rank_freale( n ) = 0;
        }
        
      }
      
    }
    
    for ( int p = 0; p < 3; p++){
      
      total += razao_rank_freale( p );
      
      //Referencia para formula do teste de Kruskal Wallis: https://pt.wikipedia.org/wiki/Teste_de_Kruskal-Wallis 
      kw( j - 1 ) = ( 12 / ( n_linha * ( n_linha + 1 ) ) ) * total - 3 * ( n_linha + 1 );
      
      //Formata cada elemento do vetor "kw" com três digitos de precisão após a virgula
      kw( j - 1 ) = roundf( kw( j - 1 ) * 1000 ) / 1000;
      
    }
    
  }
  
  return kw;
  
}
