// code.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.

#include "crypto.h"
#include <ctime>

int cmd(const char* str, int argc, char** argv) {
	for (int i = 0; i < argc; i++)
		if (strcmp(str, argv[i]) == 0) return i;
	return -1;
}

int main(int argc, char **argv) {
	srand(unsigned int(time(NULL)));
	Crypto C;
	std::list <std::string> flags{ "-v", "-de", "-kl",
		"-h", "-com", "-tk", "-f", "-ha", "-i", "-c",
		"-a", "-k", "-o", "-p", "-g"};
	int index = 0, par = 0;
	bool encrypt = true;
	bool generate = true;
	std::vector <int> args(5, -1);
	std::vector <char*> files{ "testText.txt", "cryptText.txt",
		"alphabet.raw", "vignerKey", "decrypt.txt" };
	std::list<std::string>::iterator it = flags.begin();
	for (it; it != flags.end(); it++) {
		index = cmd(it->c_str(), argc, argv);
		if (index == -1) continue;
		if (*it == "-v") par = 1;		//Vigner's alg
		else if (*it == "-h") par = 2;	//Hill's alg
		else if (*it == "-com") par = 3;	//Compare two files
		else if (*it == "-tk") args[0] = index; //Kasiski's test
		else if (*it == "-de") args[1] = 1; //Decrypt
		else if (*it == "-kl") args[2] = index; //Key lenght
		else if (*it == "-f") args[3] = index; // Frequent Analisys
		else if (*it == "-ha") args[4] = index; //Open text analisys
		else if (*it == "-i") files[0] = argv[index + 1]; //Input file file name
		else if (*it == "-c") files[1] = argv[index + 1]; //Crypt text file name
		else if (*it == "-a") files[2] = argv[index + 1]; //Alphabet file name
		else if (*it == "-k") files[3] = argv[index + 1]; //Key file file name
		else if (*it == "-o") files[4] = argv[index + 1]; //Output file file name
		else if (*it == "-g") generate = false;			//Generate key for Hill's alg or not
		else if (*it == "-p") par = 4;
		else std::cerr << "Wrong arguments\n-v to use vigner alg\n-de to decrypt" << std::endl;
	}
	if (par == 1)
		if (args[1] == -1)
			return C.encryptVigner(files[0], files[1], files[3]);
		else
			return C.decryptVigner(files[1], files[4], files[3], files[2]);
	else if (par == 2 && args[2] != 0)
		if (args[1] == -1)
			if (generate) return C.encryptHill(files[0], files[1], atoi(argv[args[2] + 1]));
			else return C.encryptHill(files[0], files[1], files[3], atoi(argv[args[2] + 1]));
		else
			return C.decryptHill(files[1], files[4], files[3], files[2], atoi(argv[args[2] + 1]));
	else if (par == 3)
		return C.compare(files[0], files[4]);
	else if (par == 4) {
		C.getText(files[0], true);
		C.printAlph(C.alphabet);
	}
	else if (args[0] != -1)
		return C.testKasiski(files[1], atoi(argv[args[0] + 1]));
	else if (args[3] != -1)
		return C.frequentAnalisys(atoi(argv[args[2] + 1]), files[0], files[1], files[3]);
	else if (args[4] != -1) {
		C.openTextAnalisys(atoi(argv[args[2] + 1]), 0, files[0], files[1], files[3]);
		C.printKey(atoi(argv[args[2] + 1]));
		return C.keyToFile(files[3]);
	}
	else
		std::cout << "Wrong argument" << std::endl;
	//return C.frequentAnalisys(5, files[0], files[1], files[3]);
	//C.encryptHill("testText.txt", "cryptText.txt", 8);
	//C.decryptHill("cryptText.txt", "decrypt.txt", "key", "alphabet.raw", 8);
	//C.decryptHill("cryptText.txt", "decrypt.txt", "key", "alphabet.raw", 8);
	//20 3 12 9 10 24 19 1 15 23 22 10 20 11 24 2 15 2 1 3 20 14 12 17 7 24 20 10 3 6 24 6 21 5 13 8 5 17 24 17 13 0 0 5 2 9 23 7 18 5 23 22 22 6 3 3 10 2 8 23 1 15 18 11 
}