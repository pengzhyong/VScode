#include "predict.h"
#include "io.h"

#include <stdio.h>

#include <string>
#include <set>
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
using namespace std;


//��Ҫ��ɵĹ��������
void predict_server(char * info[MAX_INFO_NUM], char * data[MAX_DATA_NUM], int data_num, char * filename)
{
	char* ch_flavors = info[2];
	//int flavorNums = ch_flavors[0] - '0';//�����������
	//int flavorNums = stoi(ch_flavors);
	int flavorNums = 15;

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
			if(flaIndex < flavorNums)//�ж�flavor�Ƿ���Ҫ�����ĵ�flavor��
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
	
	vector<vector<int>> seqVec;//ÿ��Ԫ����һ��flavorʱ������
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
	// ��Ҫ���������
	char * result_file = (char *)"17\n\n0 8 0 20";

	// ֱ�ӵ�������ļ��ķ��������ָ���ļ���(ps��ע���ʽ����ȷ�ԣ�����н⣬��һ��ֻ��һ�����ݣ��ڶ���Ϊ�գ������п�ʼ���Ǿ�������ݣ�����֮����һ���ո�ָ���)
	write_result(result_file, filename);
}
