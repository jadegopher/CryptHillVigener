#pragma once
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <list>
#include <vector>
#include <cmath>
#include <bitset>
#include <algorithm>
#include <valarray>
#include "bigint/BigIntegerLibrary.hh"
#define BYTE 8


class Crypto
{
	class Elem
	{
	public:
		Elem(char, int);
		char elem;
		int freq;
	};

	std::list <Elem>::iterator it;

	std::string rawText;
	std::vector <int> codeText;
	std::vector <int> cryptText;
	std::vector <int> key;

	int fileOpen(char*, std::ifstream&);
	int fileOpen(char*, std::ofstream&);
	int getNum(std::string, int&);
	int getKey(char* fileName);
	BigInteger detM(std::vector <int>, int);
	std::vector <BigInteger> adjM(std::vector <int>, int);
	std::vector <int> multM(std::vector <int> left, std::vector <int> right, int size);
	std::vector <int> multM(std::vector <int> left, std::vector <int> right, int size, int mod);
	int gcdExtended(BigInteger, BigInteger, BigInteger&, BigInteger&);
	void encryptHillAlg(std::vector <int>, std::vector <int> codeText, int);
	void encryptVignerAlg(int);
	void numText(std::string rawText, int, bool);
	int getSum(std::list <Crypto::Elem> alph);
	std::vector<char> getUnique(std::list <Crypto::Elem> alph);

public:
	std::list <Crypto::Elem> alphabet;
	std::list <Crypto::Elem> getAlph(std::string text);
	Crypto();
	int readAlph(char* fileName);
	int getText(char* fileName, bool build);
	int alphToFile(char* = "alphabet.raw");
	int cryptToFile(bool, char* = "cryptText.txt");
	int keyToFile(char* = "key");
	int encryptHill(char* = "testText.txt", char* = "cryptText.txt", int = 8);
	int encryptHill(char*, char*, char*, int = 8);
	int encryptVigner(char* = "testText.txt", char* = "cryptText.txt", int = 2);
	int encryptVigner(char*, char*, char*);

	int decryptHill(char* = "cryptText.txt", char* = "testText.txt",
		char* = "key", char* = "alphabet.raw", int = 8);

	int decryptVigner(char* = "cryptText.txt", char* = "testText.txt",
		char* = "key", char* = "alphabet.raw");

	int testKasiski(char* fileName, int nGramm);

	int frequentAnalisys(int keyLenght, char* = "testText.txt",
		char* = "cryptText.txt", char* = "possibleKey");

	int openTextAnalisys(int blockSize, int pos, char* text = "testText.txt",
		char* cryptText = "cryptText.txt", char* keyFile = "possibleHillKey");

	int compare(char* = "testText.txt", char* = "cryptText.txt");
	template<typename T> T gcd(T, T);
	template<typename T> T multiGcd(std::vector <T>);
	void printAlph(std::list <Crypto::Elem> alphabet);
	void printKey(int);
	template <typename T> void printM(std::vector <T> matrix, int blockSize, int alphSize);
	template <typename T> void printM(std::vector <T> matrix, int blockSize);
	int keyGen(int);
	~Crypto();
};

