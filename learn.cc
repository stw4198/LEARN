#include "learn.h"

int main(int argc, char** argv){

  std::string argv_str(argv[0]);
  std::string base = argv_str.substr(0, argv_str.find_last_of("/"));
  const char* path = base.c_str();

  std::vector<const char*> inputs = input_select(argc,argv);
  const char* function = inputs[0];
  const char* file = inputs[1];
  const char* component = inputs[2];
  const char* signal = inputs[3];
  const char* rPMT = inputs[4];
  const char* dTank = inputs[5];
  const char* nbins = inputs[6];
  const char* nx = inputs[7];
  const char* file_path = inputs[8];

  if(!strcmp(function,"")){
    // print statment in inputs.cc
    return -1;
  }
  if(strcmp(function,"--help")){
    printf("Inputs:\nFunction = %s\nFile = %s\nComponent = %s\nSignal = %s\nrPMT = %s mm\ndTank = %s m\nnbins = %s\nnX = %s ns\nfile path = %s\n",function,file,component,signal,rPMT,dTank,nbins,nx,file_path);
  }
  
  if(!strcmp(function,"--pdf")){
    if(!strcmp(file,"")){
      printf("No file. Please specify with -f [file]\n");
      return -1;
    }
    else if(!strcmp(component,"")){
      printf("No component. Please specify with -c [component]\n");
      return -1;
    }
    printf("\nMaking PDFs for %s\n\n\n",component);
    pdf_gen(file,component,std::stoi(nbins),1200*std::stoi(dTank),std::stoi(rPMT),nx);
  }

  else if(!strcmp(function,"--merge")){
    if(!strcmp(signal,"")){
      printf("No signal. Please specify with -s [signal]\n");
      return -1;
    }
    printf("\nMerging PDFs for %s\n\n\n",signal);
    merge_PDFs(signal);
  }

  else if(!strcmp(function,"--like")){
    if(!strcmp(file,"")){
      printf("No file. Please specify with -f [file]\n");
      return -1;
    }
    else if(!strcmp(component,"")){
      printf("No component. Please specify with -c [component]\n");
      return -1;
    }
    printf("\nCreating likelihoods for %s\n\n\n",component);
    likehood(file,component,nx);
  }

  else if(!strcmp(function,"--eval")){
    if(!strcmp(file,"")){
      printf("No file. Please specify with -f [file]\n");
      return -1;
    }
    else if(!strcmp(component,"")){
      printf("No component. Please specify with -c [component]\n");
      return -1;
    }
    const char* results_file = "results_learn.csv";
    printf("\nEvaluating likelihoods for %s\n\n\n",component);

    std::vector<double> results = likehood_classify(file,component,nx);
    double kept = results[0];
    double MC = 2*results[1];
    double rate = results[2];

    printf("\nWriting results\n\n\n");
    const char* command = Form("python3 %s/results.py %s %s %s %f %f %.9e",path,results_file,component,dTank,kept,MC,rate);
    std::system (command);
  }

  else if(!strcmp(function,"--energy")){
    if(!strcmp(signal,"")){
      printf("No signal. Please specify with -s [signal]\n");
      return -1;
    }
    printf("\nOptimising energy cuts and evaluating dwell time for %s in %s m tank\n\n\n",signal,dTank);
    const char* command = Form("python3 %s/energy_opt.py %s %s %s",path,file_path,dTank,signal);
    std::system (command);
  }

  else if(!strcmp(function,"--veto")){
    if(!strcmp(signal,"")){
      printf("No signal. Please specify with -s [signal]\n");
      return -1;
    }
    const char* results_file = "results_learn.csv";
    printf("\nEvaluating dwell time for %s in %s m tank\n\n\n",signal,dTank);
    const char* command = Form("python3 %s/veto.py %s %s %s",path,results_file,dTank,signal);
    std::system (command);
  }

  else if(!strcmp(function,"--mltrain")){
    if(!strcmp(signal,"")){
      printf("No signal. Please specify with -s [signal]\n");
      return -1;
    }
    printf("\nTraining AdaBoost for fast neutrons.\n\n\n");
    //call ml training
    const char* bg = "fn";
    const char* command = Form("python3 %s/MLEARN/fn_finder.py %s %s %s",path,file_path,signal,bg);
    std::system (command);
  }

  else if(!strcmp(function,"--mlval")){
    if(!strcmp(signal,"")){
      printf("No signal. Please specify with -s [signal]\n");
      return -1;
    }
    printf("\nApplying AdaBoost for fast neutron classification.\n\n\n");
    //call ml validation
    const char* bg = "fn";
    const char* command = Form("python3 %s/MLEARN/validation.py %s %s %s %s",path,file_path,signal,dTank,bg);
    std::system (command);
  }
  
  else if(!strcmp(function,"--help")){
    printf("\nWelcome to LEARN (Likelihood Event Analysis for Reactor Neutrinos)\n\n\n");
    printf("There are several stages to this analysis\n\n");
    printf("1) Creating PDFs\nDo learn --pdf -f [input file] -c [component]\n\n");
    printf("2) Scaling and merging PDFs\nDo learn --merge -s [signal]\n\n");
    printf("3) Creating likelihoods\nDo learn --like -f [input file] -c [component]\n\n");
    printf("4) Evaluating likelihoods and rates\nDo learn --eval -f [input file] -c [component]\n\n");
    printf("5) Data reduction\nNot currently integrated into LEARN\n\n");
    printf("6 (1)) Evaluating dwell times (energy cut)\nDo learn --energy -s [signal]\n\n");
    printf("6 (2)) Evaluating dwell times (no energy cut)\nDo learn --veto -s [signal]\n\n");
    
    std::vector<std::string> components;
    std::vector<double> rates;
    std::vector<std::string> line_values;
    std::ifstream theFile ("rates.csv");
    if(!theFile.is_open()){
      printf("File rates.csv does not exist, please copy to working directory from learn to see list of components and rates\n");
      return -1;
    }
    std::string line;
    std::getline(theFile, line);
  
    while(std::getline(theFile, line))
    {
      std::string line_value;
      std::stringstream ss(line);
      while(std::getline(ss, line_value, ','))
      {
        line_values.push_back(line_value);
      }
    }
  
    for (int i = 0; i < line_values.size(); i = i + 2) {
      components.push_back(line_values[i]);
    }
  
    for (int i = 1; i < line_values.size(); i = i + 2) {
      rates.push_back(std::stod(line_values[i]));
    }
    
    printf("Component options and rates are set in rates.csv\nThey are currently defined as:\n");
    for (int i = 0; i < components.size(); i++) {
      printf("%s: %.4e per second\n",components[i].c_str(),rates[i]);
    }
  
    
  }

  return 0;

}
