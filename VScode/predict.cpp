#include "predict.h"
//#include "io.h"
#include <stdio.h>
#include <string>
#include <set>
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <numeric>
#include <algorithm>
#define  TRAIN true
using namespace std;

void gen_seq(char * info[MAX_INFO_NUM], char * data[MAX_DATA_NUM], int data_num, vector<vector<int>>& seqVec);
void get_coef_linear(const vector<vector<int>>& seqVec, vector<vector<double>>& coefVec);//Ԥ��ϵ��
void get_coef_logistic(vector<vector<int>>& seqVec, vector<vector<double>>& coefVec);//Ԥ��ϵ��
void get_coef_linear_trend(const vector<vector<int>>& seqVec, vector<vector<double>>& coefVec, vector<vector<double>>& trendCoefVec);

void predict_seq(vector<vector<int>>& seqVec, const vector<vector<double>>& coefVec, vector<vector<int>>& predVec, int predDays);//ȫ��Ԥ����,��Ԥ������ӵ�ԭ����ĩβ
void predict_seq(vector<vector<int>>& seqVec, const vector<vector<double>>& coefVec, const vector<vector<double>>& trendCoefVec, vector<vector<int>>& predVec);

//====
void Max_CPU(int sever_num, int flavor_num[], vector<string> &sever_info);
void Max_MEM(int sever_num, int flavor_num[], vector<string> &out_info);
void Lay_flavor(int sever_num, int flavor_num, vector<vector<int>> &p, int flavor_sever[], int fla_cpu, int fla_mem);
string intostr(int a);
void InitialINT(int NEWINT[]);
void Lay_test(vector<string>& layVec, vector<string>& sever_info, bool isCPU);
int CPU_NUM = 0;
int MEM_NUM = 0;
int DISK = 0;

//��Ҫ��ɵĹ��������
void predict_server(char * info[MAX_INFO_NUM], char * data[MAX_DATA_NUM], int data_num, char * filename)
{
	string lineStr = info[0];
	istringstream line(lineStr);
	string cpuStr, memStr, diskStr, tempStr;
	line >> cpuStr >> memStr >> diskStr;

	int cpuNums = stoi(cpuStr);
	int memNums = stoi(memStr);
	int distNums = stoi(diskStr);
	CPU_NUM = cpuNums;
	MEM_NUM = memNums;
	DISK = distNums;

	lineStr = info[2];
	istringstream line1(lineStr);
	//line = istringstream(lineStr);
	line1 >> tempStr;

	int flvTypes = stoi(tempStr);//��ҪԤ���������
	vector<int> needPredictFlvs;//����������ȥ����flavor��֮�������
	for (int i = 3; i < 3 + flvTypes; i++)
	{
		lineStr = info[i];
		istringstream line(lineStr);
		//line = istringstream(lineStr);
		line >> tempStr;
		needPredictFlvs.push_back(stoi(tempStr.substr(6)));
	}
	bool isCPU = false;
	int needPredDays = 7;//��ҪԤ�������

	isCPU = (string(info[3 + flvTypes + 1]) == "CPU\n");
	string beginDate = info[3 + flvTypes + 3];
	string endDate = info[3 + flvTypes + 4];

	istringstream line2(beginDate);
	//line = istringstream(beginDate);
	line2 >> tempStr;
	tempStr = tempStr.substr(8, 2);
	int beginDay = stoi(tempStr);
	istringstream line3(endDate);
	//line = istringstream(endDate);
	line3 >> tempStr;
	tempStr = tempStr.substr(8, 2);
	int endDay = stoi(tempStr);
	needPredDays = endDay - beginDay;

	//=====ѵ��&Ԥ��
	vector<vector<int>> seqVec;//ÿ��Ԫ����һ��flavorʱ������
	gen_seq(info, data, data_num, seqVec);//ÿ��Ԫ����һ��flavorʱ������);
	vector<vector<double>> coefVec;
	vector<vector<double>> trendCoefVec;
	if (TRAIN)//��ѵ�����Ǵ��ļ���ȡѵ���õ�ϵ��
	{
		get_coef_linear(seqVec, coefVec);
		//get_coef_logistic(seqVec, coefVec);
		//get_coef_linear_trend(seqVec, coefVec, trendCoefVec);
	}
	else
	{
		ifstream inFile("coef_linear.txt");
		vector<double> oneLine;
		string line;
		while (getline(inFile, line))
		{
			oneLine.clear();
			istringstream lineStream(line);
			string word;
			while (lineStream >> word)
				oneLine.push_back(stod(word));
			coefVec.push_back(oneLine);
		}
	}
	//get_coef_logistic(seqVec, coefVec);

	//==========������
	coefVec.clear();
	for (int i = 0; i < seqVec.size(); i++)
	{
		vector<double> tempVec;
		tempVec.push_back(1);
		coefVec.push_back(tempVec);
	}

	//==========

	vector<vector<int>> predVec;
	predict_seq(seqVec, coefVec, predVec, needPredDays);
	//predict_seq(seqVec, coefVec, trendCoefVec, predVec);

	vector<int> fCnt;
	for (int i = 0; i < predVec.size(); i++)
	{
		int total = 0;
		for (auto iter = predVec.at(i).begin(); iter != predVec.at(i).end(); ++iter)
		{
			total += *iter;
		}
		if (total < 0)
			total = 0;
		fCnt.push_back(total);
	}

	//============����
	vector<int> flvNums;
	int k = 0;
	int totalFlavors = 0;
	for (auto iter = needPredictFlvs.begin(); iter != needPredictFlvs.end(); ++iter)
	{
		//flvNums.push_back(fCnt.at(*iter - 1));
		flvNums.push_back(fCnt.at(k));
		totalFlavors += fCnt.at(k);
		k++;

	}
	vector<string> layVec;
	layVec.push_back(to_string(flvTypes));
	for (int i = 0; i < flvTypes; i++)
	{
		tempStr = "flavor" + to_string(needPredictFlvs.at(i)) + " " + to_string(flvNums.at(i));
		layVec.push_back(tempStr);
	}


	vector<string> serInfoVec;
	Lay_test(layVec, serInfoVec, isCPU);

	// ��Ҫ���������
	string resultStr = to_string(totalFlavors) +"\n";
	for (auto iter = layVec.begin() + 1; iter != layVec.end(); ++iter)
	{
		resultStr += *iter + "\n";
	}
	resultStr += "\n";
	char * result_file = "";
	for (auto iter = serInfoVec.begin(); iter != serInfoVec.end(); ++iter)
	{
		istringstream line(*iter);
		//line = istringstream(*iter);
		while (line >> tempStr)
		{
			resultStr += tempStr;
			resultStr += " ";
		}
		resultStr += "\n";
	}
	resultStr.pop_back();//ȥ�����һ�����з�
	result_file = const_cast<char*>(resultStr.c_str());
	

	// ֱ�ӵ�������ļ��ķ��������ָ���ļ���(ps��ע���ʽ����ȷ�ԣ�����н⣬��һ��ֻ��һ�����ݣ��ڶ���Ϊ�գ������п�ʼ���Ǿ�������ݣ�����֮����һ���ո�ָ���)
	write_result(result_file, filename);
	system("pause");
}

void gen_seq(char * info[MAX_INFO_NUM], char * data[MAX_DATA_NUM], int data_num, vector<vector<int>>& seqVec)
{
	const char* ch_flavors = info[2];
	//int flavorNums = ch_flavors[0] - '0';//�����������
	int flavorNums = stoi(ch_flavors);
	//int flavorNums = 15;

	string* flavors = new string[flavorNums];
	for (int i = 0; i < flavorNums; i++)
	{
		flavors[i] = info[i + 3];
	}

	vector<string> flavorVec;
	for (int i = 0; i < flavorNums; i++)
	{
		istringstream sstring(flavors[i]);
		string flavorName;
		string cpus;
		string mems;
		sstring >> flavorName >> cpus >> mems;
		flavorVec.push_back(flavorName);
	}

	string beginData;
	vector<vector<int>> dataVec;//��ʷ�������ݣ�ÿһ��Ԫ��Ϊÿ�������
	dataVec.reserve(data_num);

	string iden;//���������ID
	string flavorType;//������ͺ�
	string reqDate;//��������
	int dayIndex = -1;
	string prevDate = "0000";
	for (int i = 0; i < data_num; i++)
	{
		istringstream sstring(data[i]);
		sstring >> iden >> flavorType >> reqDate;
		if (reqDate != prevDate)//������µ�һ�죬�������Ԫ��
		{
			prevDate = reqDate;
			dayIndex++;
			vector<int> dayVec(flavorNums, 0);
			dataVec.push_back(dayVec);
			int flaIndex = 0;
			string temp;
			for (auto temp : flavorVec)
			{
				if (flavorType == temp)
				{
					break;
				}
				flaIndex++;
			}
			if (flaIndex < flavorNums)//�ж�flavor�Ƿ���Ҫ�����ĵ�flavor��
				dataVec.at(dayIndex).at(flaIndex)++;
		}
		else
		{
			int flaIndex = 0;
			string temp;
			for (auto temp : flavorVec)
			{
				if (flavorType == temp)
				{
					break;
				}
				flaIndex++;
			}
			if (flaIndex < flavorNums)
				dataVec.at(dayIndex).at(flaIndex)++;
		}
	}

	//vector<vector<int>> seqVec;//ÿ��Ԫ����һ��flavorʱ������
	int dayNums = dataVec.size();
	for (int i = 0; i < flavorNums; i++)
	{
		vector<int> seqence(dayNums, 0);
		seqVec.push_back(seqence);
	}
	int daysIndex = 0;
	for (auto iter = dataVec.begin(); iter != dataVec.end(); ++iter)
	{
		for (int i = 0; i < flavorNums; i++)
		{
			seqVec.at(i).at(daysIndex) = (*iter).at(i);
		}
		daysIndex++;
	}
	//ȥ���쳣��
	for (int i = 0; i < seqVec.size(); i++)
	{
		int nums = seqVec.at(i).size();
		if (nums <= 0)
		{
			cout << "days number in negtive!" << endl;
			continue;
		}
		double average = accumulate(seqVec.at(i).begin(), seqVec.at(i).end(), 0) / (nums * 1.0);
		double variance = 0;
		for (auto iterm : seqVec.at(i))
			variance += (iterm - average)*(iterm - average);
		variance = sqrt(variance / nums);
		double posLimit = average + 3 * variance;
		double negLimit = average - 3 * variance;
		for (auto iter = seqVec.at(i).begin() + 1; iter != seqVec.at(i).end() - 1; ++iter)
		{
			if (*iter < negLimit || *iter > posLimit)
			{
				*iter = 0.5 * (*(iter - 1) + *(iter + 1));
			}
		}
	}

	//�����ݰ�����Ϊ��λ�ϲ�,Ϊ�˴������ܣ�������ǰ��ļ���
	double dropDays = dayNums % 7;
	for (auto& iterm : seqVec)
		iterm.erase(iterm.begin(), iterm.begin() + dropDays);
	vector<vector<int>> weekSeqVec;
	for (auto iter = seqVec.begin(); iter != seqVec.end(); ++iter)
	{
		vector<int> tempVec;
		int sum = 0;
		for (auto inner = iter->begin(); inner != iter->end(); inner = inner + 7)
		{
			sum = accumulate(inner, inner + 7, 0);
			tempVec.push_back(sum);
		}
		weekSeqVec.push_back(tempVec);
	}

	//����Ϊ��λ
	/*seqVec.clear();
	seqVec = weekSeqVec;*/

	string seqFileName = "ReqSeqence.txt";
	ofstream outFile(seqFileName);
	/*for (int i = 0; i < flavorNums; i++)
	{
	string temp;
	istringstream(flavors[i]) >> temp;
	outFile << temp << "\t";
	}
	outFile << endl;*/
	for (int i = 0; i < dayNums; i++)
	{
		for (int j = 0; j < flavorNums; j++)
		{
			outFile << dataVec.at(i).at(j) << "\t";
		}
		outFile << endl;
	}
	outFile.close();
}

void get_coef_linear(const vector<vector<int>>& seqVec, vector<vector<double>>& coefVec)
{
	ofstream outFile("coef_linear.txt");

	//Ԥ���㷨
	
	int coefNums = seqVec.at(0).size() > 5 ? 5 : seqVec.at(0).size();
	coefNums = 14;
	for (int flavorIndex = 0; flavorIndex < seqVec.size(); flavorIndex++)
	{
		//������ʼ��
		double* coef = new double[coefNums + 1];
		double* dw = new double[coefNums + 1];
		for (int i = 0; i < coefNums + 1; i++)
		{
			dw[i] = 0;
		}
		for (int i = 0; i < coefNums; i++)
		{
			coef[i] = 0.1 * (rand() % 10 - 5);
		} 
		coef[coefNums] = 0;
		
		//����������
		double learnRate = 0.001;
		double eps = 1;
		double delta = 1000;
		//��ʼѵ��,�ݶ��½�������ģ��
		int t = 0;
		double preDelta = 0;
		int iterCnt = 0;
		while (abs(delta) > eps)// && abs(delta - preDelta) > 0.01)
		{
			for (int i = 0; i < coefNums + 1; i++)
			{
				dw[i] = 0;
			}
			//if (iterCnt % 1000 == 0)
			//	learnRate *= 0.1;

			if (iterCnt++ > 5000)
				break;//������ѭ��
			preDelta = delta;
			delta = 0;
			for (int i = 0; i < seqVec.at(flavorIndex).size() - coefNums; i++)
			{
				//delta += seqVec.at(flavorIndex).at(i + coefNums);
				double y = seqVec.at(flavorIndex).at(i + coefNums);
				double yHat = 0;
				for (int j = 0; j < coefNums; j++)
				{
					//delta = delta - coef[j] * seqVec.at(flavorIndex).at(i + j);
					yHat += coef[j] * seqVec.at(flavorIndex).at(i + j);
				}
				//delta -= coef[coefNums];
				yHat += coef[coefNums];
				if (yHat < 0)
					yHat = 0;
				delta += 0.5 * (y - yHat) * (y - yHat);
				//lost = 0.5*(y - yHat)^2
				for (int j = 0; j < coefNums; j++)
				{
					if(yHat >=0)
						dw[j] += -(y - yHat)*seqVec.at(flavorIndex).at(i + j);

				}
				if (yHat >= 0)
					dw[coefNums] += yHat - y;
			}
			for (int i = 0; i < coefNums + 1; i++)
			{
				dw[i] = dw[i] / (seqVec.at(flavorIndex).size() - coefNums);
			}

			for (int i = 0; i < coefNums + 1; i++)
			{
				coef[i] = coef[i] - learnRate * dw[i];
			}
			cout << "the delta of " << t++ << " interators: "
				<< abs(delta) <<  "\tflavor: " << flavorIndex <<endl;
		}
		//��ӡϵ��
		for (int i = 0; i < coefNums; i++)
		{
			cout << coef[i] << ends;
			outFile << coef[i] << " ";//���ϵ��
		}
		cout << endl;

		outFile << endl;

		vector<double> tempVec;
		for (int i = 0; i < coefNums + 1; i++)
		{
			tempVec.push_back(coef[i]);
		}
		coefVec.push_back(tempVec);
		
		delete[] coef;
		delete[] dw;
	}

	outFile.close();
}

void get_coef_logistic(vector<vector<int>>& seqVec, vector<vector<double>>& coefVec)
{
	//�������ݹ�һ��,�������ݲ����й�һ�����޸�logstic��ʽ��ʹ�������ݷ�Χƥ��
	vector<vector<double>> seqVecDoub;
	vector<pair<double, double>> resumScale;//��������ϵ����bias������һ��
	for (auto it1 = seqVec.begin(); it1 != seqVec.end(); ++it1)
	{
		vector<double> tempVec;
		for (auto it2 = (*it1).begin(); it2 != (*it1).end(); ++it2)
		{
			tempVec.push_back(*it2);
		}
		seqVecDoub.push_back(tempVec);
	}

	//Ԥ���㷨
	//vector<vector<int>> seqVec;
	int coefNums = 14;
	for (int flavorIndex = 0; flavorIndex < seqVecDoub.size(); flavorIndex++)
	{
		//���ݹ�һ��
		int maxValue = *max_element(seqVecDoub.at(flavorIndex).begin(), seqVecDoub.at(flavorIndex).end());
		int minValue = *min_element(seqVecDoub.at(flavorIndex).begin(), seqVecDoub.at(flavorIndex).end());
		double scales = maxValue - minValue;
		double bias = accumulate(seqVecDoub.at(flavorIndex).begin(), seqVecDoub.at(flavorIndex).end(), 0) / (seqVecDoub.at(flavorIndex).size() * 1.0);
		//for (auto& iterm : seqVecDoub.at(flavorIndex))
			//iterm = iterm / scales + bias;

		pair<double, double> tempPair(scales, bias);
		resumScale.push_back(tempPair);//���ڷ���һ��

		//������ʼ��
		double* coef = new double[coefNums + 1];
		double* dw = new double[coefNums + 1];
		for (int i = 0; i < coefNums + 1; i++)
		{
			dw[i] = 0;
		}
		for (int i = 0; i < coefNums; i++)
		{
			coef[i] = 0.1 * (rand() % 10 - 5);
		}
		coef[coefNums] = 0;//ƫ��ϵ��b

		//����������
		double learnRate = 0.001;
		double eps = 1;
		double delta = 1000;
		//��ʼѵ��,�ݶ��½�������ģ��
		int t = 0;
		double preDelta = 0;
		int iterCnt = 0;
		while (abs(delta) > eps)// && abs(delta - preDelta) > 0.01)
		{
			for (int i = 0; i < coefNums + 1; i++)
			{
				dw[i] = 0;
			}
			if (iterCnt++ > 2000)
				break;//������ѭ��
			preDelta = delta;
			delta = 0;
			for (int i = 0; i < seqVecDoub.at(flavorIndex).size() - coefNums; i++)
			{
				double y = seqVecDoub.at(flavorIndex).at(i + coefNums);
				double z = 0;
				double yHat = 0;
				for (int j = 0; j < coefNums; j++)
				{
					z += coef[j] * seqVecDoub.at(flavorIndex).at(i + j);
				}
				z += coef[coefNums]; 
				//z /= scales;
				//yHat = scales / (1 + exp(-z)) + bias;
				//yHat = 1.0 / (1.0 + exp(-z));

				//ReLU
				yHat = z > 0 ? z : 0;

				delta += 0.5 * (y - yHat) * (y - yHat);
				//lost = -ylogyHat - (1 - y)log(1 - yHat)
				for (int j = 0; j < coefNums; j++)
				{
					if(z > 0)
						dw[j] += (yHat - y)*seqVecDoub.at(flavorIndex).at(i + j);

				}
				dw[coefNums] += yHat - y;
			}
			for (int i = 0; i < coefNums + 1; i++)
			{
				dw[i] = dw[i] / (seqVecDoub.at(flavorIndex).size() - coefNums);
				//if (abs(dw[i]) > 3)
				//	dw[i] = 0;
			}
			
			for (int i = 0; i < coefNums + 1; i++)
			{
				coef[i] = coef[i] - learnRate * dw[i];
			}
			cout << "the delta of " << t++ << " interators: "
				<< abs(delta) <<"\tflavors: "<< flavorIndex + 1 << endl;
		}
		//��ӡϵ��
		/*for (int i = 0; i < coefNums; i++)
		{
			cout << coef[i] << ends;
		}
		cout << endl;*/

		vector<double> tempVec;
		for (int i = 0; i < coefNums + 1; i++)
		{
			tempVec.push_back(coef[i]);
		}
		coefVec.push_back(tempVec);

		delete[] coef;
		delete[] dw;
	}


}

void get_coef_linear_trend(const vector<vector<int>>& seqVec, vector<vector<double>>& coefVec, vector<vector<double>>& trendCoefVec)//����������
{
	ofstream outFile("coef_linear.txt");
	ofstream trendOutFile("coef_linear_trend.txt");

	//Ԥ���㷨
	int coefNums = seqVec.at(0).size() > 5 ? 5 : seqVec.at(0).size();
	coefNums = 5;
	for (int flavorIndex = 0; flavorIndex < seqVec.size(); flavorIndex++)
	{
		//������ʼ��
		double* coef = new double[coefNums + 1];
		double* dw = new double[coefNums + 1];
		double* trendCoef = new double[coefNums];
		double* trendDw = new double[coefNums];
		for (int i = 0; i < coefNums + 1; i++)
		{
			dw[i] = 0;
		}
		for (int i = 0; i < coefNums; i++)
		{
			coef[i] = 0.1 * (rand() % 10 - 5);
			trendDw[i] = 0;
			trendCoef[i] = 0.1 * (rand() % 10 - 5);
		}
		coef[coefNums] = 0;

		//����������
		double learnRate = 0.0001;
		double eps = 1;
		double delta = 1000;
		//��ʼѵ��,�ݶ��½�������ģ��
		int t = 0;
		double preDelta = 0;
		int iterCnt = 0;
		while (abs(delta) > eps)// && abs(delta - preDelta) > 0.01)
		{

			for (int i = 0; i < coefNums + 1; i++)
			{
				dw[i] = 0;
			}
			if (iterCnt++ > 2000)
				break;//������ѭ��
			preDelta = delta;
			delta = 0;
			for (int i = 0; i < seqVec.at(flavorIndex).size() - coefNums; i++)
			{
				//delta += seqVec.at(flavorIndex).at(i + coefNums);
				double y = seqVec.at(flavorIndex).at(i + coefNums);
				double yHat = 0;
				for (int j = 0; j < coefNums; j++)
				{
					//delta = delta - coef[j] * seqVec.at(flavorIndex).at(i + j);
					int preSeq = seqVec.at(flavorIndex).at(i + j);
					if(i + j > 0)
						preSeq = seqVec.at(flavorIndex).at(i + j - 1);
					yHat += coef[j] * seqVec.at(flavorIndex).at(i + j) + trendCoef[j] * (seqVec.at(flavorIndex).at(i + j) - preSeq);
				}
				//delta -= coef[coefNums];
				yHat += coef[coefNums];
				//if (yHat < 0)
				//	yHat = 0;
				delta += 0.5 * (y - yHat) * (y - yHat);
				//lost = 0.5*(y - yHat)^2
				
				//����dw
				for (int j = 0; j < coefNums; j++)
				{
					//if (yHat >= 0)
					dw[j] += -(y - yHat)*seqVec.at(flavorIndex).at(i + j);
					int preSeq = seqVec.at(flavorIndex).at(i + j);
					if (i + j > 0)
						preSeq = seqVec.at(flavorIndex).at(i + j - 1);
					trendDw[i] += -(y - yHat)*(seqVec.at(flavorIndex).at(i + j) - preSeq);
				}
				//if (yHat >= 0)
				dw[coefNums] += yHat - y;
			}
			for (int i = 0; i < coefNums + 1; i++)
			{
				dw[i] = dw[i] / (seqVec.at(flavorIndex).size() - coefNums);
				if( i < coefNums)
					trendDw[i] = trendDw[i] / (seqVec.at(flavorIndex).size() - coefNums);
			}

			for (int i = 0; i < coefNums + 1; i++)
			{
				coef[i] = coef[i] - learnRate * dw[i];
				if (i < coefNums)
					trendCoef[i] -= learnRate * trendDw[i];
			}
			cout << "the delta of " << t++ << " interators: "
				<< abs(delta) << "\tflavor: " << flavorIndex << endl;
		}
		//��ӡϵ��
		for (int i = 0; i < coefNums; i++)
		{
			cout << coef[i] << ends;
			outFile << coef[i] << " ";//���ϵ��
			trendOutFile << trendCoef[i] << " ";
		}
		cout << endl;

		outFile << endl;
		trendOutFile << endl;

		vector<double> tempVec;
		vector<double> trendVec;
		for (int i = 0; i < coefNums + 1; i++)
		{
			tempVec.push_back(coef[i]);
			if (i < coefNums)
				trendVec.push_back(trendCoef[i]);
		}
		coefVec.push_back(tempVec);
		trendCoefVec.push_back(trendVec);

		delete[] coef;
		delete[] dw;
		delete[] trendDw;
		delete[] trendCoef;
	}

	outFile.close();
}

void predict_seq(vector<vector<int>>& seqVec, const vector<vector<double>>& coefVec, vector<vector<int>>& predVec, int predDays)
{
		
	//int predNums = seqVec.at(0).size() > 1 ? 1 : seqVec.at(0).size();
	int predNums = predDays;
	int coefNums = coefVec.at(0).size() - 1;
	int seqLen = seqVec.at(0).size();
	for (int i = 0; i < coefVec.size(); i++)
	{
		vector<int> tempVec;
		for (int j = 0; j < predNums; j++)
		{
			double predValue = 0;
			for (int k = 0; k < coefNums; k++)
			{
				//double c = coefVec.at(i).at(k);
				//double v = seqVec.at(i).at(seqLen - 1 - coefNums + k + j);
				predValue += coefVec.at(i).at(k) * (seqVec.at(i).at(seqLen - coefNums + k + j));
			}
			predValue += coefVec.at(i).at(coefNums);

			seqVec.at(i).push_back(static_cast<int>(predValue + 0.5));
			tempVec.push_back(static_cast<int>(predValue + 0.5));
		}
		predVec.push_back(tempVec);
	}
}

void predict_seq(vector<vector<int>>& seqVec, const vector<vector<double>>& coefVec, const vector<vector<double>>& trendCoefVec, vector<vector<int>>& predVec)
{

	//int predNums = seqVec.at(0).size() > 1 ? 1 : seqVec.at(0).size();
	int predNums = 1;
	int coefNums = coefVec.at(0).size() - 1;
	int seqLen = seqVec.at(0).size();
	for (int i = 0; i < coefVec.size(); i++)
	{
		vector<int> tempVec;
		for (int j = 0; j < predNums; j++)
		{
			double predValue = 0;
			for (int k = 0; k < coefNums; k++)
			{
				//double c = coefVec.at(i).at(k);
				//double v = seqVec.at(i).at(seqLen - 1 - coefNums + k + j);
				double preSeq = seqVec.at(i).at(seqLen - coefNums + k + j);
				if (seqLen - coefNums + k + j > 0)
					preSeq = seqVec.at(i).at(seqLen - coefNums + k + j - 1);
				predValue += coefVec.at(i).at(k) * (seqVec.at(i).at(seqLen - coefNums + k + j))
					+ trendCoefVec.at(i).at(k) * (seqVec.at(i).at(seqLen - coefNums + k + j) - preSeq);//������
			}
			predValue += coefVec.at(i).at(coefNums);

			seqVec.at(i).push_back(static_cast<int>(predValue + 0.5));
			tempVec.push_back(static_cast<int>(predValue + 0.5));
		}
		predVec.push_back(tempVec);
	}
}

//========
void Lay_test(vector<string>& layVec, vector<string>& sever_info,bool isCPU)
{
	string flavor1 = "flavor1";
	string flavor2 = "flavor2";
	string flavor3 = "flavor3";
	string flavor4 = "flavor4";
	string flavor5 = "flavor5";
	string flavor6 = "flavor6";
	string flavor7 = "flavor7";
	string flavor8 = "flavor8";
	string flavor9 = "flavor9";
	string flavor10 = "flavor10";
	string flavor11 = "flavor11";
	string flavor12 = "flavor12";
	string flavor13 = "flavor13";
	string flavor14 = "flavor14";
	string flavor15 = "flavor15";

	int Sever_num = 0;//����ķ���������

	int i;
	int flavor_num[15] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };//����flavor����
	int cpu_num = 0;//CPU��������
	int mem_num = 0;//�ڴ���������(GB)
	for (i = 1; i < layVec.size(); i++)
	{
		string line = layVec.at(i);
		int line_length = line.length();
		int first_blank = line.find_first_of(" ");
		int last_blank = line.find_last_of(" ");
		string VM_name = line.substr(0, first_blank);
		int VM_num = atoi(line.substr(last_blank).c_str());
		if (!VM_name.compare(flavor1))
		{
			flavor_num[0] = VM_num;
		}
		else if (!VM_name.compare(flavor2))
		{
			flavor_num[1] = VM_num;
		}
		else if (!VM_name.compare(flavor3))
		{
			flavor_num[2] = VM_num;
		}
		else if (!VM_name.compare(flavor4))
		{
			flavor_num[3] = VM_num;
		}
		else if (!VM_name.compare(flavor5))
		{
			flavor_num[4] = VM_num;
		}
		else if (!VM_name.compare(flavor6))
		{
			flavor_num[5] = VM_num;
		}
		else if (!VM_name.compare(flavor7))
		{
			flavor_num[6] = VM_num;
		}
		else if (!VM_name.compare(flavor8))
		{
			flavor_num[7] = VM_num;
		}
		else if (!VM_name.compare(flavor9))
		{
			flavor_num[8] = VM_num;
		}
		else if (!VM_name.compare(flavor10))
		{
			flavor_num[9] = VM_num;
		}
		else if (!VM_name.compare(flavor11))
		{
			flavor_num[10] = VM_num;
		}
		else if (!VM_name.compare(flavor12))
		{
			flavor_num[11] = VM_num;
		}
		else if (!VM_name.compare(flavor13))
		{
			flavor_num[12] = VM_num;
		}
		else if (!VM_name.compare(flavor14))
		{
			flavor_num[13] = VM_num;
		}
		else if (!VM_name.compare(flavor15))
		{
			flavor_num[14] = VM_num;
		}
	}

	cpu_num = flavor_num[0] * 1 + flavor_num[1] * 1 + flavor_num[2] * 1 + flavor_num[3] * 2 + flavor_num[4] * 2
		+ flavor_num[5] * 2 + flavor_num[6] * 4 + flavor_num[7] * 4 + flavor_num[8] * 4 + flavor_num[9] * 8
		+ flavor_num[10] * 8 + flavor_num[11] * 8 + flavor_num[12] * 16 + flavor_num[13] * 16 + flavor_num[14] * 16;

	mem_num = flavor_num[0] * 1 + flavor_num[1] * 2 + flavor_num[2] * 4 + flavor_num[3] * 2 + flavor_num[4] * 4
		+ flavor_num[5] * 8 + flavor_num[6] * 4 + flavor_num[7] * 8 + flavor_num[8] * 16 + flavor_num[9] * 8
		+ flavor_num[10] * 16 + flavor_num[11] * 32 + flavor_num[12] * 16 + flavor_num[13] * 32 + flavor_num[14] * 64;

	Sever_num = max(cpu_num / CPU_NUM, mem_num / MEM_NUM) + 1;//�������������һ��Ԥ�����������

	if (isCPU)
	{
		Max_CPU(Sever_num, flavor_num, sever_info);
	}
	else
	{
		Max_MEM(Sever_num, flavor_num, sever_info);
	}
	for (i = 0; i < sever_info.size(); i++)
	{
		cout << sever_info[i] << endl;
	}
}

void Max_CPU(int sever_num, int flavor_num[], vector<string> &out_info)
{
	vector<vector<int>> p;//����������
	int i, j;
	for (i = 0; i < sever_num; i++)
	{
		vector<int> sever_x = { i + 1,CPU_NUM,MEM_NUM };//������������{ID��CPU,�ڴ�}
		p.push_back(sever_x);//�������ܼ���
	}

	//����flavor15
	//for (i = 0; i < flavor_num[14]; i++)
	//{
	//	int flavor15_sever[1000] = {};//������flavor15�ķ�������Ϣ
	//	for (j=0;j<sever_num;j++)
	//	{
	//		//���CPU���ڴ��㹻�������
	//		if (p[j][1] > 16 && p[j][2] > 64)
	//		{
	//			p[j][1] -= 16;//��j+1��������CPU���ڴ����
	//			p[j][2] -= 64;
	//			flavor15_sever[j] += 1;//������j+1���������ϵ�flavor15����+1
	//			break;
	//		}
	//		//������һ��������Ҳ�治�£����½�һ��������������
	//		else if (j==sever_num-1)
	//		{
	//			vector<int> sever_info = { sever_num + 1,CPU_NUM,MEM_NUM };//������sever_num+1����������{ID��CPU,�ڴ�}
	//			p.push_back(sever_info);//��ӵ��������ܼ���
	//			sever_num++;//����������+1
	//			p[j][1] -= 16;
	//			p[j][2] -= 64;
	//			flavor15_sever[j] += 1;
	//			break;
	//		}
	//	}
	//}
	int flavor15_sever[1000];//������flavor15�ķ���������Щ
	InitialINT(flavor15_sever);
	Lay_flavor(sever_num, flavor_num[14], p, flavor15_sever, 16, 64);

	int flavor14_sever[1000];
	InitialINT(flavor14_sever);
	Lay_flavor(sever_num, flavor_num[13], p, flavor14_sever, 16, 32);

	int flavor13_sever[1000];
	InitialINT(flavor13_sever);
	Lay_flavor(sever_num, flavor_num[12], p, flavor13_sever, 16, 16);

	int flavor12_sever[1000];
	InitialINT(flavor12_sever);
	Lay_flavor(sever_num, flavor_num[11], p, flavor12_sever, 8, 32);

	int flavor11_sever[1000];
	InitialINT(flavor11_sever);
	Lay_flavor(sever_num, flavor_num[10], p, flavor11_sever, 8, 16);

	int flavor10_sever[1000];
	InitialINT(flavor10_sever);
	Lay_flavor(sever_num, flavor_num[9], p, flavor10_sever, 8, 8);

	int flavor9_sever[1000];
	InitialINT(flavor9_sever);
	Lay_flavor(sever_num, flavor_num[8], p, flavor9_sever, 4, 16);

	int flavor8_sever[1000];
	InitialINT(flavor8_sever);
	Lay_flavor(sever_num, flavor_num[7], p, flavor8_sever, 4, 8);

	int flavor7_sever[1000];
	InitialINT(flavor7_sever);
	Lay_flavor(sever_num, flavor_num[6], p, flavor7_sever, 4, 4);

	int flavor6_sever[1000];
	InitialINT(flavor6_sever);
	Lay_flavor(sever_num, flavor_num[5], p, flavor6_sever, 2, 8);

	int flavor5_sever[1000];
	InitialINT(flavor5_sever);
	Lay_flavor(sever_num, flavor_num[4], p, flavor5_sever, 2, 4);

	int flavor4_sever[1000];
	InitialINT(flavor4_sever);
	Lay_flavor(sever_num, flavor_num[3], p, flavor4_sever, 2, 2);

	int flavor3_sever[1000];
	InitialINT(flavor3_sever);
	Lay_flavor(sever_num, flavor_num[2], p, flavor3_sever, 1, 4);

	int flavor2_sever[1000];
	InitialINT(flavor2_sever);
	Lay_flavor(sever_num, flavor_num[1], p, flavor2_sever, 1, 2);

	int flavor1_sever[1000];
	InitialINT(flavor1_sever);
	Lay_flavor(sever_num, flavor_num[0], p, flavor1_sever, 1, 1);

	out_info.push_back(intostr(sever_num));
	for (i = 0; i < sever_num; i++)
	{
		string info = "";
		info = intostr(i + 1);
		//����÷���������flavor15
		if (flavor15_sever[i] > 0)
		{
			info += " flavor15 " + intostr(flavor15_sever[i]);
		}
		if (flavor14_sever[i] > 0)
		{
			info += " flavor14 " + intostr(flavor14_sever[i]);
		}
		if (flavor13_sever[i] > 0)
		{
			info += " flavor13 " + intostr(flavor13_sever[i]);
		}
		if (flavor12_sever[i] > 0)
		{
			info += " flavor12 " + intostr(flavor12_sever[i]);
		}
		if (flavor11_sever[i] > 0)
		{
			info += " flavor11 " + intostr(flavor11_sever[i]);
		}
		if (flavor10_sever[i] > 0)
		{
			info += " flavor10 " + intostr(flavor10_sever[i]);
		}
		if (flavor9_sever[i] > 0)
		{
			info += " flavor9 " + intostr(flavor9_sever[i]);
		}
		if (flavor8_sever[i] > 0)
		{
			info += " flavor8 " + intostr(flavor8_sever[i]);
		}
		if (flavor7_sever[i] > 0)
		{
			info += " flavor7 " + intostr(flavor7_sever[i]);
		}
		if (flavor6_sever[i] > 0)
		{
			info += " flavor6 " + intostr(flavor6_sever[i]);
		}
		if (flavor5_sever[i] > 0)
		{
			info += " flavor5 " + intostr(flavor5_sever[i]);
		}
		if (flavor4_sever[i] > 0)
		{
			info += " flavor4 " + intostr(flavor4_sever[i]);
		}
		if (flavor3_sever[i] > 0)
		{
			info += " flavor3 " + intostr(flavor3_sever[i]);
		}
		if (flavor2_sever[i] > 0)
		{
			info += " flavor2 " + intostr(flavor2_sever[i]);
		}
		if (flavor1_sever[i] > 0)
		{
			info += " flavor1 " + intostr(flavor1_sever[i]);
		}
		out_info.push_back(info);
	}
}

void Max_MEM(int sever_num, int flavor_num[], vector<string> &out_info)
{
	vector<vector<int>> p;//����������
	int i, j;
	for (i = 0; i < sever_num; i++)
	{
		vector<int> sever_x = { i + 1,CPU_NUM,MEM_NUM };//������������{ID��CPU,�ڴ�}
		p.push_back(sever_x);//�������ܼ���
	}

	int flavor15_sever[1000];//������flavor15�ķ���������Щ
	InitialINT(flavor15_sever);
	Lay_flavor(sever_num, flavor_num[14], p, flavor15_sever, 16, 64);

	int flavor14_sever[1000];
	InitialINT(flavor14_sever);
	Lay_flavor(sever_num, flavor_num[13], p, flavor14_sever, 16, 32);

	int flavor12_sever[1000];
	InitialINT(flavor12_sever);
	Lay_flavor(sever_num, flavor_num[11], p, flavor12_sever, 8, 32);

	int flavor13_sever[1000];
	InitialINT(flavor13_sever);
	Lay_flavor(sever_num, flavor_num[12], p, flavor13_sever, 16, 16);


	int flavor11_sever[1000];
	InitialINT(flavor11_sever);
	Lay_flavor(sever_num, flavor_num[10], p, flavor11_sever, 8, 16);

	int flavor9_sever[1000];
	InitialINT(flavor9_sever);
	Lay_flavor(sever_num, flavor_num[8], p, flavor9_sever, 4, 16);

	int flavor10_sever[1000];
	InitialINT(flavor10_sever);
	Lay_flavor(sever_num, flavor_num[9], p, flavor10_sever, 8, 8);

	int flavor8_sever[1000];
	InitialINT(flavor8_sever);
	Lay_flavor(sever_num, flavor_num[7], p, flavor8_sever, 4, 8);

	int flavor6_sever[1000];
	InitialINT(flavor6_sever);
	Lay_flavor(sever_num, flavor_num[5], p, flavor6_sever, 2, 8);

	int flavor7_sever[1000];
	InitialINT(flavor7_sever);
	Lay_flavor(sever_num, flavor_num[6], p, flavor7_sever, 4, 4);


	int flavor5_sever[1000];
	InitialINT(flavor5_sever);
	Lay_flavor(sever_num, flavor_num[4], p, flavor5_sever, 2, 4);

	int flavor3_sever[1000];
	InitialINT(flavor3_sever);
	Lay_flavor(sever_num, flavor_num[2], p, flavor3_sever, 1, 4);

	int flavor4_sever[1000];
	InitialINT(flavor4_sever);
	Lay_flavor(sever_num, flavor_num[3], p, flavor4_sever, 2, 2);

	int flavor2_sever[1000];
	InitialINT(flavor2_sever);
	Lay_flavor(sever_num, flavor_num[1], p, flavor2_sever, 1, 2);

	int flavor1_sever[1000];
	InitialINT(flavor1_sever);
	Lay_flavor(sever_num, flavor_num[0], p, flavor1_sever, 1, 1);

	out_info.push_back(intostr(sever_num));
	for (i = 0; i < sever_num; i++)
	{
		string info = "";
		info = intostr(i + 1);
		//����÷���������flavor15
		if (flavor15_sever[i] > 0)
		{
			info += " flavor15 " + intostr(flavor15_sever[i]);
		}
		if (flavor14_sever[i] > 0)
		{
			info += " flavor14 " + intostr(flavor14_sever[i]);
		}
		if (flavor13_sever[i] > 0)
		{
			info += " flavor13 " + intostr(flavor13_sever[i]);
		}
		if (flavor12_sever[i] > 0)
		{
			info += " flavor12 " + intostr(flavor12_sever[i]);
		}
		if (flavor11_sever[i] > 0)
		{
			info += " flavor11 " + intostr(flavor11_sever[i]);
		}
		if (flavor10_sever[i] > 0)
		{
			info += " flavor10 " + intostr(flavor10_sever[i]);
		}
		if (flavor9_sever[i] > 0)
		{
			info += " flavor9 " + intostr(flavor9_sever[i]);
		}
		if (flavor8_sever[i] > 0)
		{
			info += " flavor8 " + intostr(flavor8_sever[i]);
		}
		if (flavor7_sever[i] > 0)
		{
			info += " flavor7 " + intostr(flavor7_sever[i]);
		}
		if (flavor6_sever[i] > 0)
		{
			info += " flavor6 " + intostr(flavor6_sever[i]);
		}
		if (flavor5_sever[i] > 0)
		{
			info += " flavor5 " + intostr(flavor5_sever[i]);
		}
		if (flavor4_sever[i] > 0)
		{
			info += " flavor4 " + intostr(flavor4_sever[i]);
		}
		if (flavor3_sever[i] > 0)
		{
			info += " flavor3 " + intostr(flavor3_sever[i]);
		}
		if (flavor2_sever[i] > 0)
		{
			info += " flavor2 " + intostr(flavor2_sever[i]);
		}
		if (flavor1_sever[i] > 0)
		{
			info += " flavor1 " + intostr(flavor1_sever[i]);
		}
		out_info.push_back(info);
	}
}
//���������
void Lay_flavor(int sever_num, int flavor_num, vector<vector<int>> &p, int flavor_sever[], int fla_cpu, int fla_mem)
{
	int i, j;
	//��������flavor
	for (i = 0; i < flavor_num; i++)
	{

		for (j = 0; j < sever_num; j++)
		{
			//���CPU���ڴ��㹻�������
			if (p[j][1] > fla_cpu && p[j][2] > fla_mem)
			{
				p[j][1] -= fla_cpu;//��j+1��������CPU���ڴ����
				p[j][2] -= fla_mem;
				flavor_sever[j] += 1;//������j+1���������ϵ�flavorX����+1
				break;
			}
			//������һ��������Ҳ�治�£����½�һ��������������
			else if (j == sever_num - 1)
			{
				vector<int> sever_info = { sever_num + 1,CPU_NUM,MEM_NUM };//������sever_num+1����������{ID��CPU,�ڴ�}
				p.push_back(sever_info);//��ӵ��������ܼ���
				sever_num++;//����������+1
				p[j][1] -= fla_cpu;
				p[j][2] -= fla_mem;
				flavor_sever[j] += 1;//������j+1���������ϵ�flavorX����+1
				break;
			}
		}

	}
}

//intת��string����
string intostr(int a)
{
	ostringstream stream;
	stream << a;  //nΪint����
	return stream.str();
}

//�������ʼ��Ϊ0������1000
void InitialINT(int NEWINT[])
{
	int i;
	for (i = 0; i < 1000; i++)
	{
		NEWINT[i] = 0;
	}
}