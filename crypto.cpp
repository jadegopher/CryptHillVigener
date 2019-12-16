#include "crypto.h"

Crypto::Crypto() {

}

Crypto::Elem::Elem(char elem, int freq) {
	Elem::elem = elem;
	Elem::freq = freq;
}

//Find in string number and convert to int from current position
//Returns number and saves position at first symbol after a number
int Crypto::getNum(std::string str, int &pos) {
	std::string number;
	for (pos; std::isdigit(str[pos]); pos++)
		number += str[pos];
	//std::cout << number << std::endl;
	return std::stoi(number);
}

//Open file as binary and returns size of file. Returns -1 if file doesn't exist
int Crypto::fileOpen(char* fileName, std::ifstream &inFile) {
	inFile.open(fileName, std::ios::binary | std::ios::ate);
	if (!inFile.is_open()) {
		std::cerr << "File does not find";
		return -1;
	}
	int size = int(inFile.tellg());
	inFile.seekg(0, std::ios::beg);
	//std::cout << "Total file lenght: " << size << std::endl;
	return size;
}

int Crypto::fileOpen(char* fileName, std::ofstream& outFile) {
	outFile.open(fileName, std::ios::binary);
	if (!outFile.is_open()) {
		std::cerr << "File does not find";
		return -1;
	}
	return 1;
}

std::list <Crypto::Elem> Crypto::getAlph(std::string text) {
	std::list <Crypto::Elem> alph;
	for (size_t i = 0; i < text.size(); i++) {
		for (it = alph.begin(); it != alph.end(); it++)
			if (it->elem == text[i]) {
				it->freq++;
				break;
			}
		if (it == alph.end()) {
			//std::cout << "New element: \"" << rawFile[i] << "\" has appended. Code:" <<
			//	+unsigned char(rawFile[i]) << std::endl;
			alph.push_back(Elem(text[i], 1));
		}
	}
	alph.sort([](const Elem& a, const Elem& b) {return a.freq > b.freq; });
	return alph;
}

//Gets alphabet from text file 
//and calculate frequnce of repeating each letter
//Returns -1 if file does not found
int Crypto::getText(char* fileName, bool build) {
	std::ifstream inFile;
	//Get text from file
	int size = fileOpen(fileName, inFile);
	if (size == -1) return -1;
	rawText = std::string((std::istreambuf_iterator<char>(inFile)),
		(std::istreambuf_iterator<char>()));
	inFile.close();
	if (build) {
		alphabet.clear();
		alphabet = getAlph(rawText);
	}
	return 1;
}

//Read key from File
int Crypto::getKey(char* fileName) {
	std::ifstream inFile;
	int size = fileOpen(fileName, inFile);
	if (size == -1) return -1;
	key = std::vector <int>();
	for (int j = 0; inFile >> j;) key.push_back(j);
	if (key.size() == 0) return -1;
	inFile.close();
	return 1;
}

//Stores alphaped to direct file.
//Stores to alphabet.raw by default
int Crypto::alphToFile(char* fileName) {
	std::ofstream outFile;
	if (fileOpen(fileName, outFile) == -1) return -1;
	for (it = alphabet.begin(); it != alphabet.end(); it++)
		outFile << it->elem << " " << it->freq << " ";
	outFile.close();
	return 1;
}

//Build text by using cryptText vector and
//Stores crypted text into a file
//If flag en is false determinate from last 8 bytes which part of text 
//should to delete (for Hill's alg)
int Crypto::cryptToFile(bool en, char* fileName) {
	std::ofstream outFile;
	if (fileOpen(fileName, outFile) == -1) return -1;
	int del = 0;
	if (!en) {
		std::string tmp;
		for (int i = cryptText.size() - BYTE; i < int(cryptText.size()); i++)
			tmp += std::to_string(cryptText[i]);
		//std::cout << tmp << std::endl;
		del = std::stoi(tmp, 0, 2);
	}
	for (int i = 0, j = 0; i < int(cryptText.size()) - del; i++, j = 0)
		for (it = alphabet.begin(); it != alphabet.end(); it++, j++)
			if (cryptText[i] == j) {
				outFile << it->elem;
				break;
			}
	outFile.close();
	return 1;
}

//writes key to file
int Crypto::keyToFile(char* fileName) {
	std::ofstream outFile;
	if (fileOpen(fileName, outFile) == -1) return -1;
	for (int i = 0; i < int(key.size()); i++)
		outFile << key[i] << " ";
	outFile.close();
	return 1;
}

//Get alphabet from file (and frequencies)
//Stores to alphabet list
int Crypto::readAlph(char* fileName) {
	std::ifstream inFile;
	if(fileOpen(fileName, inFile) == -1)
		return -1;
	std::string rawAlph((std::istreambuf_iterator<char>(inFile)),
		(std::istreambuf_iterator<char>()));
	inFile.close();
	alphabet.clear();
	char tmp;
	for (int i = 0; i < int(rawAlph.size()); i++) {
		tmp = rawAlph[i];
		alphabet.push_back(Elem(tmp, getNum(rawAlph, i += 2)));
	}
	return 1;
}

//Hill's encrypt algorithm
void Crypto::encryptHillAlg(std::vector <int> lKey, std::vector <int> codeText, int blockSize) {
	cryptText = std::vector <int>(codeText.size());
	int powMod = blockSize;
	if(blockSize != 1)
	powMod = int(pow(blockSize, 2)) - 1;
	//std::cout << "CRYPT: " << codeText.size() << std::endl;
	for (int i = 0, j = 0, sum = 0, blockStart = 0, k = 0;
		true;
		i++, j += blockSize) {
		if (j >= powMod + 1) {
			j %= powMod;
			cryptText[k] = sum % alphabet.size();
			//std::cout << cryptText[k] << " ";
			//std::cout << std::endl;
			if (j == blockSize) {
				blockStart += blockSize;
				if (blockStart == codeText.size()) break;
				j = 0;
			}
			sum = 0;
			k++;
			i = blockStart;
		}
		//std::cout << codeText[i] << " * " << lKey[j] << std::endl;
		sum += codeText[i] * lKey[j];
		//std::cout << i << " ";
	}
	//std::cout << std::endl;
}

void Crypto::encryptVignerAlg(int blockSize) {
	cryptText = std::vector <int>(codeText.size());
	for (int i = 0, j = 0; i < int(codeText.size()); i++, j++) {
		if (j == blockSize) j = 0;
		cryptText[i] = (codeText[i] + key[j]) % alphabet.size();
	}
}

//Get source text. Build an alphabet. Crypt text and write to result file
int Crypto::encryptHill(char *source, char *result, int blockSize) {
	if (getText(source, true) == -1) return -1;
	if (alphToFile() == -1) return -2;
	numText(rawText, blockSize, true);
	if (keyGen(blockSize) == -1) return -3;
	printKey(blockSize);
	encryptHillAlg(key, codeText, blockSize);
	if (cryptToFile(true, result) == -1) return -4;
	if (keyToFile() == -1) return -5;
	return 1;
}

//Do the same as cryptHill but you could specify key
int Crypto::encryptHill(char* source, char* result, char* keyFile, int blockSize) {
	if (getText(source, true) == -1) return -1;
	if (alphToFile() == -1) return -2;
	numText(rawText, blockSize, true);
	if (getKey(keyFile) == -1) return -3;
	if (key.size() % blockSize != 0) {
		std::cerr << "Bad key";
		key = std::vector <int>();
		return -1;
	}
	printKey(blockSize);
	encryptHillAlg(key, codeText, blockSize);
	if (cryptToFile(true, result) == -1) return -4;
	return 1;
}

//Open file. Generates alphabet. Numerate all symbols
//and encrypt text
int Crypto::encryptVigner(char* source, char* result, int keySize) {
	if (getText(source, true) == -1) return -1;
	if (alphToFile() == -1) return -2;
	numText(rawText, 0, false);
	if (keyGen(keySize) == -1) return -3;
	printKey(keySize);
	encryptVignerAlg(keySize);
	if (cryptToFile(true, result) == -1) return -4;
	if (keyToFile() == -1) return -5;
	return 1;
}

//DO the same but u can upload key from direct file
int Crypto::encryptVigner(char* source, char* result, char* keyFile) {
	if (getText(source, true) == -1) return -1;
	if (alphToFile() == -1) return -2;
	numText(rawText, 0, false);
	if (getKey(keyFile) == -1) return -3;
	encryptVignerAlg(key.size());
	if (cryptToFile(true, result) == -1) return -4;
	return 1;
}

int Crypto::decryptHill(char* source, char* result, char* keyFile, char* alphFile, int blockSize) {
	if (blockSize < 1) return -4;
	if (getText(source, true) == -1) return -1;
	if (readAlph(alphFile) == -1) return -6;
	numText(rawText, blockSize, false);
	if (getKey(keyFile) == -1) return -3;
	//printAlph(alphabet);

	BigInteger det = detM(key, blockSize);
	//std::cout << "Determinate is: " << det << std::endl;
	BigInteger x, y;
	int aSize = alphabet.size();
	int g = gcdExtended(det, aSize, x, y);
	//std::cout << "GCD: " << g << "\nX: " << x << "\nY: " << y << std::endl;
	if (g == 0) return -5;
	if (det > 0 && x < 0) x += BigInteger(aSize);
	//std::cout << "GCD: " << g << "\nX: " << x << "\nY: " << y << std::endl;
	//printKey(blockSize);
	std::vector <BigInteger> bAdj = adjM(key, blockSize);
	std::vector <int> modAdj(key.size());
	BigInteger tmp, temp;
	for (int i = 0, j = 0; i < int(bAdj.size()); i++, j += blockSize) {
		if(j >= int(bAdj.size())) j %= bAdj.size() - 1;

		tmp = bAdj[i] % aSize;
		temp = (tmp * x) % aSize;
		//std::cout << "DETERMINANT: " << bAdj[i] << " DMOD: " << 
		//	tmp << " DMODX: " << tmp * x << " DMODXMOD: " << temp << std::endl;
		modAdj[j] = temp.toInt();
	}
	//printM(modAdj, blockSize);
	//printM(multM(modAdj, key, blockSize), blockSize, aSize);
	encryptHillAlg(modAdj, codeText, blockSize);
	if (cryptToFile(false, result) == -1) return -7;
	return 1;
}

//Vigner's decription
//Notice that cryptText is decrypt text!!!
int Crypto::decryptVigner(char* source, char* result, char* keyFile, char* alphFile) {
	if (getText(source, true) == -1) return -1;
	if (readAlph(alphFile) == -1) return -2;
	numText(rawText, 0, false);
	if (getKey(keyFile) == -1) return -3;
	printKey(key.size());
	cryptText = std::vector <int>(codeText.size());
	for (int i = 0, j = 0; i < int(codeText.size()); i++, j++) {
		if (j == key.size()) j = 0;
		cryptText[i] = (codeText[i] + alphabet.size() - key[j]) % alphabet.size();
	}
	if (cryptToFile(true, result) == -1) return -4;
	return 1;
}

//Test Kasiski 
//Returns the most freequent interval
int Crypto::testKasiski(char* fileName, int nGramm) {
	if (getText(fileName, true) == -1) return -1;
	if (int(rawText.size()) < nGramm) return -2; //Source text is too small
	std::string tmp(rawText.substr(0, nGramm));
	std::vector <int> interval;
	std::map <int, int> gcd;
	for (int pos = 0, find = 0; pos < int(rawText.size());) {
		find = rawText.find(tmp, find + 1);
		if (find != std::string::npos)
			//std::cout << find << " " << pos << std::endl;
			interval.push_back(find - pos);
		else {
			if (interval.size() != 0) {
				auto g = multiGcd(interval);
				auto it = gcd.find(g);
				if (it != gcd.end())
					it->second++;
				else
					gcd.insert(std::pair <int,int>(g, 1));
				interval = std::vector <int>();
			}
			pos++;
			find = pos;
			tmp = std::string(rawText.substr(pos, nGramm));
		}
	}
	std::pair <int, int> max(1, 1);
	for (auto i = gcd.begin(); i != gcd.end(); i++)
		//std::cout << i->first << " " << i->second << std::endl;
		if (i->second > max.second) {
			max.first = i->first;
			max.second = i->second;
		}
	std::cout << "Max value is: " << max.first << " " << max.second << std::endl;
	return max.first;
}

int Crypto::openTextAnalisys(int blockSize, int pos, char* text,
	char* cryptText, char* keyFile) {
	if (getText(text, true) == -1) return -1;
	numText(rawText, blockSize, false);
	std::vector<int> textSample;
	BigInteger det, x, y;
	int g = 0, aSize = alphabet.size(), mod = 0;
	size_t in = pos;
	//printAlph(alphabet);
	textSample = std::vector <int>();
	for (size_t j = in; j < codeText.size() && j < in + pow(blockSize, 2); j++)
		textSample.push_back(codeText[j]);
	det = detM(textSample, blockSize);
	g = gcdExtended(det, aSize, x, y);
	y = (det * x) % aSize;
	if (g == 0) return -5;
	x = (x % aSize + aSize) % aSize;
	
	std::vector <BigInteger> bAdj = adjM(textSample, blockSize);
	std::vector <int> modAdj(textSample.size());
	BigInteger tmp, temp;
	for (int i = 0, j = 0; i < int(bAdj.size()); i++, j += blockSize) {
		if (j >= int(bAdj.size())) j %= bAdj.size() - 1;

		tmp = bAdj[i] % aSize;
		temp = (tmp * x) % aSize;
		modAdj[j] = temp.toInt();
	}
	if (getText(cryptText, false) == -1) return -2;
	numText(rawText, blockSize, false);
	std::vector <int> cryptTextSample;

	for (size_t i = in; i < in + pow(blockSize, 2); i++) {
		cryptTextSample.push_back(codeText[i]);
	}

	key = std::vector <int>(int(pow(blockSize, 2)));
	key = multM(cryptTextSample, modAdj, blockSize, aSize);
	//printM(multM(key, textSample, blockSize, aSize), blockSize);
	//printM(cryptTextSample, blockSize);
	if (multM(key, textSample, blockSize, aSize) == cryptTextSample) {
		det = detM(key, blockSize);
		g = gcdExtended(det, aSize, x, y);
		y = (det * x) % alphabet.size();
		mod = y.toInt();
		if (mod != 1) 
			openTextAnalisys(blockSize, in += size_t(pow(blockSize, 2)), text, cryptText, keyFile);
	}
	else 
		openTextAnalisys(blockSize, in += size_t(pow(blockSize, 2)), text, cryptText, keyFile);
	return 1;
}


std::vector<char> Crypto::getUnique(std::list <Crypto::Elem> alph) {
	std::map <int, int> unique;
	std::vector <char> ret;
	for (auto i = alph.begin(); i != alph.end(); i++) {
		auto it = unique.find(i->freq);
		if (it != unique.end())
			it->second++;
		else
			unique.insert(std::pair<int, int>(i->freq, 1));
	}
	for (auto i = alph.begin(); i != alph.end(); i++)
		for (auto it = unique.begin(); it != unique.end(); it++)
			if (i->freq == it->first && it->second == 1) ret.push_back(i->elem);
	return ret;
}

//Gets alphabet from file
//Parse crypt text by blocks e.t. 
//0 5 10 15...
//1 6 11 16...
//2 7 12 17... if key lenght 5
//and using frequent analisys calculates key
//and stores into key variable
//and saves in file
int Crypto::frequentAnalisys(int keyLenght,
	char* testText, char* cryptTextName, char* keyFile) {
	if (getText(cryptTextName, false) == -1) return -1;
	std::list <Crypto::Elem> cryptAlph;
	std::string cryptText = rawText;
	cryptAlph = getAlph(cryptText);
	if (getText(testText, true) == -1) return -2;
	std::list <Crypto::Elem> miniAlph;
	std::list <Crypto::Elem> miniOpenAlph;
	key = std::vector <int>();
	int k = 0;
	int l = 0;
	for (int i = 0; i < keyLenght; i++) {
		std::string miniText;
		std::string miniOpenText;
		for (size_t j = i; j < cryptText.size(); j += keyLenght)
			miniText += cryptText[j];
		for (size_t j = i; j < rawText.size(); j += keyLenght)
			miniOpenText += rawText[j];

		miniAlph = getAlph(miniText);
		miniOpenAlph = getAlph(miniOpenText);
		std::vector <char> unique = getUnique(miniOpenAlph);
		std::vector <char> cryptUnique = getUnique(miniAlph);
		if (cryptUnique.size() > 0){
			k = 0;
			l = 0;
			int res = -1, tmp = -1;
			for (auto it = alphabet.begin(); it != alphabet.end(); it++, k++, l++) {
				if (it->elem == unique[0])
					res = l;
				if (it->elem == cryptUnique[0])
					tmp = k;
				if (res != -1 && tmp != -1) {
					key.push_back(abs(tmp - res));
					break;
				}
			}
		}
		else {
			k = 0;
			for (auto it = alphabet.begin(); it != alphabet.end(); it++, k++) {
				if (it->elem == miniAlph.begin()->elem) {
					key.push_back(k);
					break;
				}
			}
		}
	}
	printKey(keyLenght);
	keyToFile(keyFile);
	return 1;
}

//Gets sum of all symbols in alphabet
int Crypto::getSum(std::list <Crypto::Elem> alph) {
	int ret = 0;
	for (auto it = alph.begin(); it != alph.end(); it++)
		ret += it->freq;
	return ret;
}

//Numerates all symbols from source file
//Saves result in codeText vector
//If true, adds to end some bytes. The last 8 bytes
//shows how many bytes we had added (for Hill's alg)
void Crypto::numText(std::string rawText, int blockSize, bool en) {
	//printAlph(alphabet);
	int rawTextSize = rawText.size();
	codeText = std::vector <int>(rawTextSize);
	int size = rawTextSize;
	if (en) {
		if (BYTE % blockSize != 0) size += BYTE - (BYTE % blockSize) + blockSize;
		else size += BYTE;
		codeText = std::vector <int>(size);
		if (size % blockSize != 0) {
			size = size - (size % blockSize) + blockSize;
			codeText = std::vector <int>(size);
		}
	}
	//std::cout << "\nNUM: " << rawTextSize << " " << blockSize << " " << size << std::endl;
	for (int i = 0, j = 0; i < rawTextSize; i++, j = 0)
		for (it = alphabet.begin(); it != alphabet.end(); it++, j++)
			if (rawText[i] == it->elem) {
				codeText[i] = j;
				//std::cout << codeText[i] << " ";
				break;
			}
	if (en) {
		for (int i = rawTextSize; i < size - BYTE; i++) {
			codeText[i] = 0;
			//std::cout << codeText[i] << " ";
		}
		std::bitset <BYTE> delBin = size - rawTextSize;
		for (int i = size - BYTE, j = BYTE - 1; i < size; i++, j--) {
			codeText[i] = delBin[j];
			//std::cout << codeText[i] << " ";
		}
		//std::cout << std::endl;
	}
}

//Generates key and saves it in key vector
int Crypto::keyGen(int blockSize) {
	if (alphabet.size() <= 1 || blockSize < 1) return -1;
	key = std::vector <int> (int(pow(blockSize, 2)));
	for (int i = 0; i < int(key.size()); i++) 
		key[i] = rand() % alphabet.size();
	BigInteger det = detM(key, blockSize);
	//std::cout << "Determinate is: " << det << std::endl;
	BigInteger x, y;
	int aSize = alphabet.size();
	int g = gcdExtended(det, aSize, x, y);
	y = (det * x) % alphabet.size();
	int mod = y.toInt();
	//std::cout << mod << std::endl;
	if (mod != 1) keyGen(blockSize);
	/*for (int i = 0; i < int(codeText.size()); i++) {
		if (i % blockSize == 0) std::cout << std::endl;
		std::cout << codeText[i] << " ";
	}
	std::cout << "\n" << std::endl;*/
	return 1;
}

//Calculates determinant of a matrix
BigInteger Crypto::detM(std::vector <int> inM, int n) {
	int index; // Initialize result
	BigInteger det(1), total(1), num1, num2;
	if (n == 0) return det;
	// temporary array for storing row  
	std::vector <BigInteger> matrix(inM.size());
	std::vector <BigInteger> temp(n);
	for (int i = 0; i < int(matrix.size()); i++)
		matrix[i] = inM[i];
	//loop for traversing the diagonal elements  
	for (int i = 0; i < n; i++) {
		index = i; // initialize the index   
		//finding the index which has non zero value   
		while (index < n && matrix[index * n + i] == 0) index++;
		if (index == n) continue;
		if (index != i) {
			//loop for swaping the diagonal element row and index row   
			for (int j = 0; j < n; j++)
				std::swap(matrix[index * n + j], matrix[i * n + j]);
			det = det * int(pow(-1, index - i));
		}
		//storing the values of diagonal row elements   
		for (int j = 0; j < n; j++) temp[j] = matrix[i * n + j];
		//traversing every row below the diagonal element   
		for (int j = i + 1; j < n; j++)
		{
			num1 = temp[i]; //value of diagonal element   
			num2 = matrix[j * n + i]; //value of next row element   
			//traversing every column of row   
			// and multiplying to every row   
			for (int k = 0; k < n; k++)	
				matrix[j * n + k] = (num1 * matrix[j * n + k]) - (num2 * temp[k]);
			total = total * num1; // Det(kA)=kDet(A);   
		}
	}
	//mulitplying the diagonal elements to get determinant   
	for (int i = 0; i < n; i++)	det = det * matrix[i * n + i];
	return (det / total);  
}

//Finds matrix of alg adjections
std::vector <BigInteger> Crypto::adjM(std::vector <int> inM, int n) {
	if (n == 1) return std::vector <BigInteger>(inM[0]);
	std::vector <int> matrix(int(pow(n - 1, 2)));
	std::vector <BigInteger> adj(inM.size());
	//std::cout << matrix.size() << std::endl;
	for (int i = 0, j = 0, row = 0, col = 0, rowM = 0, colM = 0; true ;j++) {
		if (j == n) {
			j = 0;
			i++;
		}
		if (i == n) {
			i = 0, rowM = 0, colM = 0, j = 0;
			//std::cout << "Det for elem: " 
			//	<< inM[row * n + col] << " " << detM(matrix, n - 1) << std::endl;
			adj[row * n + col] = BigInteger(int(pow(-1, row + col))) * detM(matrix, n - 1);
			//std::cout << adj[row * n + col] << " ";
			col++;
			if (col == n) {
				col = 0;
				row++;
				//std::cout << std::endl;
				if (row == n) break;
			}
		}
		if (i != row && j != col) {
			matrix[rowM * (n - 1) + colM] = inM[i * n + j];
			colM++;
			if (colM == n - 1) {
				colM = 0;
				rowM++;
			}
		}
	}
	return adj;
}

//Multiplicates two matricies
std::vector <int> Crypto::multM(std::vector <int> left,
	std::vector <int> right, int size) {
	std::vector <int> res;
	int tmp = 0;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			for (int k = 0; k < size; k++)
				tmp += right[i * size + k] * left[k * size + j];
			res.push_back(tmp);
			tmp = 0;
		}
	}
	return res;
}

//Multiplicates two matricies
std::vector <int> Crypto::multM(std::vector <int> left,
	std::vector <int> right, int size, int mod) {
	std::vector <int> res;
	int tmp = 0;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			for (int k = 0; k < size; k++)
				tmp += right[i * size + k] * left[k * size + j];
			res.push_back(tmp % mod);
			tmp = 0;
		}
	}
	return res;
}


// Function for extended Euclidean Algorithm  
int Crypto::gcdExtended(BigInteger a, BigInteger b,
	BigInteger &x, BigInteger &y) {
	if (a == 0){
		x = 0;
		y = 1;
		return b.toInt();
	}
	x = a;
	y = b;
	BigInteger tmp1, tmp2, temp;
	BigInteger arr[4] = { 1, 0, 0, 1 };
	while (y != 0) {
		tmp1 = x / y;
		tmp2 = x % y;
		x = y;
		y = tmp2;
		temp = arr[0];
		arr[0] = arr[1];
		arr[1] = temp - arr[1] * tmp1;
		temp = arr[2];
		arr[2] = arr[3];
		arr[3] = temp - arr[3] * tmp1;
	}
	a = x;
	x = arr[0];
	y = arr[2];
	return a.toInt();
}

template <typename T> T Crypto::gcd(T a, T b) {
	T r;
	while (b != 0) {
		r = a % b;
		a = b;
		b = r;
	}
	return a;
}

template <typename T> T Crypto::multiGcd(std::vector <T> arr) {
	int result = arr[0];
	for (int i = 1; i < arr.size(); i++)
		result = gcd(arr[i], result);
	return result;
}

int Crypto::compare(char* file, char* another) {
	std::ifstream file1, file2;
	if (fileOpen(file, file1) == -1) return -1;
	if (fileOpen(another, file2) == -1) return -2;
	std::string tmp((std::istreambuf_iterator<char>(file1)),
		(std::istreambuf_iterator<char>()));
	std::string tmp1((std::istreambuf_iterator<char>(file2)),
		(std::istreambuf_iterator<char>()));
	file1.close();
	file2.close();
	int cmp = tmp.compare(tmp1);
	if (cmp == 0) {
		std::cout << "OK" << std::endl;
		return cmp;
	}
	else
		return cmp;
}


//Just prints alphabet
void Crypto::printAlph(std::list <Crypto::Elem> alphabet) {
	std::cout << "Symbol" << "\t" << "Symbol's code" << 
		"\t" << "Frequency" << "\t" << "Number" << std::endl;
	int i = 0;
	for (it = alphabet.begin(); it != alphabet.end(); it++, i++)
		std::cout << it->elem << "\t" << +it->elem << "\t" << 
		it->freq << "\t" << i << std::endl;
}

//Just prints key
void Crypto::printKey(int blockSize) {
	std::cout << "Key: " ;
	for (int i = 0; i < int(key.size()); i++) {
		if (i % blockSize == 0) std::cout << std::endl;
		std::cout << key[i] << " ";
	}
	std::cout << std::endl;
}

//Just prints matrix
template <typename T>
void Crypto::printM(std::vector <T> matrix, int blockSize, int alphSize) {
	std::cout << "Matrix: ";
	for (int i = 0; i < int(matrix.size()); i++) {
		if (i % blockSize == 0) std::cout << std::endl;
		std::cout << matrix[i] % alphSize << " ";
	}
	std::cout << std::endl;
}

//Just prints matrix
template <typename T>
void Crypto::printM(std::vector <T> matrix, int blockSize) {
	std::cout << "Matrix: ";
	for (int i = 0; i < int(matrix.size()); i++) {
		if (i % blockSize == 0) std::cout << std::endl;
		std::cout << matrix[i] << " ";
	}
	std::cout << std::endl;
}

Crypto::~Crypto() {

}