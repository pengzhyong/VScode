#include "predict.h"
#include "io.h"

#include <stdio.h>

#include <string>
#include <set>
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <numeric>
using namespace std;

void gen_seq(char * info[MAX_INFO_NUM], char * data[MAX_DATA_NUM], int data_num, vector<vector<int>>& seqVec);
void get_coef_linear(const vector<vector<int>>& seqVec, vector<vector<double>>& coefVec);//Ԥ��ϵ��
void get_coef_logistic(const vector<vector<int>>& seqVec, vector<vector<double>>& coefVec);//Ԥ��ϵ��

void predict_seq(vector<vector<int>>& seqVec, const vector<vector<double>>& coefVec, vector<vector<int>>& predVec);//ȫ��Ԥ����,��Ԥ������ӵ�ԭ����ĩβ

//��Ҫ��ɵĹ��������
void predict_server(char * info[MAX_INFO_NUM], char * data[MAX_DATA_NUM], int data_num, char * filename)
{
	vector<vector<int>> seqVec;//ÿ��Ԫ����һ��flavorʱ������
	gen_seq(info, data, data_num, seqVec);//ÿ��Ԫ����һ��flavorʱ������);
	vector<vector<double>> coefVec;
	get_coef_linear(seqVec, coefVec);
	//get_coef_logistic(seqVec, coefVec);

	vector<vector<int>> predVec;
	predict_seq(seqVec, coefVec, predVec);

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

	// ��Ҫ���������
	char * result_file = (char *)"17\n\n0 8 0 20";

	// ֱ�ӵ�������ļ��ķ��������ָ���ļ���(ps��ע���ʽ����ȷ�ԣ�����н⣬��һ��ֻ��һ�����ݣ��ڶ���Ϊ�գ������п�ʼ���Ǿ�������ݣ�����֮����һ���ո�ָ���)
	write_result(result_file, filename);
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
	//Ԥ���㷨
	//vector<vector<int>> seqVec;
	
	int coefNums = 7;
	for (int flavorIndex = 0; flavorIndex < seqVec.size(); flavorIndex++)
	{
		//������ʼ��
		double* coef = new double[coefNums];
		double* dw = new double[coefNums];
		for (int i = 0; i < coefNums; i++)
		{
			dw[i] = 0;
		}
		for (int i = 0; i < coefNums; i++)
		{
			coef[i] = 0.1 * (rand() % 10 - 5);
		} 

		
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
			for (int i = 0; i < coefNums; i++)
			{
				dw[i] = 0;
			}
			if (iterCnt++ > 1000)
				break;//������ѭ��
			preDelta = delta;
			delta = 0;
			for (int i = 0; i < seqVec.at(flavorIndex).size() - coefNums; i++)
			{
				delta += seqVec.at(flavorIndex).at(i + coefNums);
				double y = seqVec.at(flavorIndex).at(i + coefNums);
				double yHat = 0;
				for (int j = 0; j < coefNums; j++)
				{
					delta = delta - coef[j] * seqVec.at(flavorIndex).at(i + j);
					yHat += coef[j] * seqVec.at(flavorIndex).at(i + j);
				}
				//lost = 0.5*(y - yHat)^2
				//delta += (y - yHat) * (-)
				for (int j = 0; j < coefNums; j++)
				{
					dw[j] += -(y - yHat)*seqVec.at(flavorIndex).at(i + j);

				}
			}
			for (int i = 0; i < coefNums; i++)
			{
				dw[i] = dw[i] / (seqVec.at(flavorIndex).size() - coefNums);
			}

			int kk = 1;
			for (int i = 0; i < coefNums; i++)
			{
				coef[i] = coef[i]
					- learnRate * dw[i];// accumulate(seqVec.at(flavorIndex).begin() + coefNums + i, seqVec.at(flavorIndex).end() - coefNums + i, 0) * posNeg;
			}
			cout << "the delta of " << t++ << " interators: "
				<< abs(delta) << endl;
		}
		//��ӡϵ��
		for (int i = 0; i < coefNums; i++)
		{
			cout << coef[i] << ends;
		}
		cout << endl;

		vector<double> tempVec;
		for (int i = 0; i < coefNums; i++)
		{
			tempVec.push_back(coef[i]);
		}
		coefVec.push_back(tempVec);
		
		delete[] coef;
		delete[] dw;
	}

	
}

void get_coef_logistic(const vector<vector<int>>& seqVec, vector<vector<double>>& coefVec)
{
	//Ԥ���㷨
	//vector<vector<int>> seqVec;

	int coefNums = 14;
	for (int flavorIndex = 0; flavorIndex < seqVec.size(); flavorIndex++)
	{
		//������ʼ��
		double* coef = new double[coefNums];
		double* dw = new double[coefNums];
		for (int i = 0; i < coefNums; i++)
		{
			dw[i] = 0;
		}
		for (int i = 0; i < coefNums; i++)
		{
			coef[i] = 0.001 * (rand() % 10 - 5);
		}


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
			for (int i = 0; i < coefNums; i++)
			{
				dw[i] = 0;
			}
			if (iterCnt++ > 2000)
				break;//������ѭ��
			preDelta = delta;
			delta = 0;
			for (int i = 0; i < seqVec.at(flavorIndex).size() - coefNums; i++)
			{
				delta += seqVec.at(flavorIndex).at(i + coefNums);
				double y = seqVec.at(flavorIndex).at(i + coefNums);
				int z = 0;
				double yHat = 0;
				for (int j = 0; j < coefNums; j++)
				{
					delta = delta - coef[j] * seqVec.at(flavorIndex).at(i + j);
					z += coef[j] * seqVec.at(flavorIndex).at(i + j);
				}
				yHat = 1 / (1 + exp(-z)) - 0.5;
				//lost = 0.5*(y - yHat)^2
				for (int j = 0; j < coefNums; j++)
				{
					dw[j] += (yHat - y)*seqVec.at(flavorIndex).at(i + j);

				}
			}
			for (int i = 0; i < coefNums; i++)
			{
				dw[i] = dw[i] / (seqVec.at(flavorIndex).size() - coefNums);
				if (abs(dw[i]) > 3)
					dw[i] = 0;
			}

			int kk = 1;
			for (int i = 0; i < coefNums; i++)
			{
				coef[i] = coef[i]
					- learnRate * dw[i];// accumulate(seqVec.at(flavorIndex).begin() + coefNums + i, seqVec.at(flavorIndex).end() - coefNums + i, 0) * posNeg;
				if (abs(coef[i]) > 5)
				{
					coef[i] = 0.1 * (rand() % 10) - 5;
				}
			}
			cout << "the delta of " << t++ << " interators: "
				<< abs(delta) << endl;
		}
		//��ӡϵ��
		for (int i = 0; i < coefNums; i++)
		{
			cout << coef[i] << ends;
		}
		cout << endl;

		vector<double> tempVec;
		for (int i = 0; i < coefNums; i++)
		{
			tempVec.push_back(coef[i]);
		}
		coefVec.push_back(tempVec);

		delete[] coef;
		delete[] dw;
	}


}

void predict_seq(vector<vector<int>>& seqVec, const vector<vector<double>>& coefVec, vector<vector<int>>& predVec)
{
	int predNums = 14;
	int coefNums = coefVec.at(0).size();
	int seqLen = seqVec.at(0).size();
	for (int i = 0; i < coefVec.size(); i++)
	{
		vector<int> tempVec;
		for (int j = 0; j < predNums; j++)
		{
			double predValue = 0;
			for (int k = 0; k < coefNums; k++)
			{
				double c = coefVec.at(i).at(k);
				double v = seqVec.at(i).at(seqLen - 1 - predNums + k + j);
				predValue += coefVec.at(i).at(k) * (seqVec.at(i).at(seqLen - 1 - predNums + k + j));
			}
			seqVec.at(i).push_back(static_cast<int>(predValue + 0.5));
			tempVec.push_back(static_cast<int>(predValue + 0.5));
		}
		predVec.push_back(tempVec);
	}
}
