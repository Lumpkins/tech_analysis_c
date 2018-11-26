#include "include.h"
#include <math.h>


double max(PyObject * list){
	int size;
	size=PyList_GET_SIZE(list);
	double maximum;
	

	for (int i=0;i<size;i++){
		double current=PyFloat_AsDouble(PyList_GetItem(list,i));
		if (i==0){
			maximum=current;
		}
		if ((current>maximum)|| i==0){
			maximum=current;
		}
	}

	return maximum;
}

double min(PyObject * list){
	int size;
	size=PyList_GET_SIZE(list);
	double minimum;
	

	for (int i=0;i<size;i++){
		double current=PyFloat_AsDouble(PyList_GetItem(list,i));
		if (i==0){
			minimum=current;
		}
		if ((current<minimum) || (i==0)){
			minimum=current;
		}
	}

	return minimum;
}

double sum(PyObject * list){
	int size;
  	size=PyList_GET_SIZE(list);

	double sum=0;
	
	for (int i=0;i<size;i++){
		sum=sum+PyFloat_AsDouble(PyList_GetItem(list,i));

	}

	return sum;
}

PyObject* money_flow_index(PyObject* self, PyObject* args)
{
	PyObject * TP;

	

	PyObject * RMF;

	int period=14;

	if (!PyArg_ParseTuple(args,"O!O!",&PyList_Type,&TP,&PyList_Type,&RMF))
	{
		goto error;
	}

    int size;
	size=PyList_GET_SIZE(TP);

	PyObject *lst = PyList_New(size);
	//allocate the array
	for (int i=0;i<period;i++){
		
		PyObject *num = PyFloat_FromDouble(0);
		PyList_SET_ITEM(lst, i, num); 
	}

	for (int i=period;i<size;i++){
		double pos=0;
		double neg=0;	
		for (int k = i-period+1;k<i;k++){		
			double misc=PyFloat_AsDouble(PyList_GetItem(TP,k));
			double misc2=PyFloat_AsDouble(PyList_GetItem(TP,k-1));
		
			if (misc>misc2){
				pos=pos+ PyFloat_AsDouble( PyList_GetItem(RMF,k));
			}
			else
			{
				neg=neg+ PyFloat_AsDouble( PyList_GetItem(RMF,k));
			}	
		}
		PyObject *num = PyFloat_FromDouble(100-100/(1+pos/neg));
		PyList_SET_ITEM(lst, i, num); 		
	}

	return lst;
error:
	return 0;
}

PyObject* william_r(PyObject* self, PyObject* args)
{
	PyObject * high;
	PyObject * low;
	PyObject * close;

	int period=14;

	if (!PyArg_ParseTuple(args,"O!O!O!",&PyList_Type,&high,&PyList_Type,&low,&PyList_Type,&close))
	{
		goto error;
	}//relocate the specification context

    int size;	
	size=PyList_GET_SIZE(high);

	PyObject *R = PyList_New(size);
	//allocate the array
	for (int i=0;i<period;i++){
		PyObject *num = PyFloat_FromDouble(0);
		PyList_SET_ITEM(R, i, num); 
	}

	for (int i=period;i<size;i++){
		PyObject * high_m ;
		high_m= PyList_GetSlice(high,i-period,i);
		PyObject * low_m ;
		low_m= PyList_GetSlice(low,i-period,i);

		double misc;
		double close_current;

		close_current=PyFloat_AsDouble(PyList_GetItem(close,i));

		misc=-100.*(max(high_m)-close_current)/(max(high_m)-min(low_m));

		PyObject *num = PyFloat_FromDouble(misc);
		PyList_SET_ITEM(R, i, num); 		
	}

	return R;
error:
	return 0;
}

PyObject* ultimate_ocillator(PyObject* self, PyObject* args)
{
	PyObject * high;

	PyObject * low;

	PyObject * close;

	int period=28;

	if (!PyArg_ParseTuple(args,"O!O!O!",&PyList_Type,&high,&PyList_Type,&low,&PyList_Type,&close))
	{
		goto error;
	}

    int size;
	
	size=PyList_Size(high);

	PyObject *UO = PyList_New(size);

	PyObject *BP = PyList_New(size);
	PyObject *TR = PyList_New(size);

	PyList_SetItem(UO, 0, PyFloat_FromDouble(0)); 
	PyList_SetItem(BP, 0, PyFloat_FromDouble(0)); 
	PyList_SetItem(TR, 0, PyFloat_FromDouble(0)); 

	for (int i=1;i<size;i++){
		double low_current;
		low_current=PyFloat_AsDouble(PyList_GetItem(low,i));

		double high_current;
		high_current=PyFloat_AsDouble(PyList_GetItem(high,i));

		double close_current;
		close_current=PyFloat_AsDouble(PyList_GetItem(close,i));

		double close_prev;
		close_prev=PyFloat_AsDouble(PyList_GetItem(close,i-1));

		double minimum;
		if (low_current<close_prev){
			minimum=low_current;
		}
		else{
			minimum=close_prev;
		}

		double maximum;
		if (high_current>close_prev){
			maximum=high_current;
		}
		else{
			maximum=close_prev;
		}

		PyList_SetItem(BP, i, PyFloat_FromDouble(close_current-minimum)); 
		PyList_SetItem(TR, i, PyFloat_FromDouble(maximum-minimum)); 

		double Avg7=0 ;
		double Avg14 =0;
		double Avg28 =0;

		if (i>=7){
			double s1=sum(PyList_GetSlice(BP,i-7,i));
			double s2=sum(PyList_GetSlice(TR,i-7,i)); 

			Avg7=s1/s2;
		}
		if (i>=14){
			double s1=sum(PyList_GetSlice(BP,i-14,i));
			double s2=sum(PyList_GetSlice(TR,i-14,i));

			Avg14=s1/s2;
		}
		if (i>=28){
			double s1=sum(PyList_GetSlice(BP,i-28,i));
			double s2=sum(PyList_GetSlice(TR,i-28,i));

			Avg28=s1/s2;
		}

		PyObject *num = PyFloat_FromDouble(100*(4*Avg7+2*Avg14+Avg28)/7);

		PyList_SetItem(UO, i, num); 		
	}

	return  UO;
error:
	return 0;
}

PyObject* find_spans_one_above_other(PyObject* self, PyObject* args)
{
	PyObject * s1;

	PyObject * s2;

	PyObject * index;
	
	if (!PyArg_ParseTuple(args,"O!O!O!",&PyList_Type,&s1,&PyList_Type,&s2,&PyList_Type,&index))
	{
		goto error;
	}

    int size;
	
	size=PyList_GET_SIZE(s1);

	PyObject *pairs ;
		
	pairs=PyList_New(0);

	PyObject * item;

	int i_save;
	i_save=0;

	bool s1gts2;
	
	//https://docs.python.org/3.6/extending/extending.html#a-simple-example

	if (PyFloat_AsDouble(PyList_GetItem(s1,0))>PyFloat_AsDouble(PyList_GetItem(s2,0)))
	{		
		s1gts2=true;
	}
	else
	{
		s1gts2=false;
	}
	
	for (int i=1;i<size;i++)
	{
		double s1_current;
		double s2_current;
		s1_current=PyFloat_AsDouble(PyList_GetItem(s1,i));
		s2_current=PyFloat_AsDouble(PyList_GetItem(s2,i));
	
		if((s1_current>s2_current) && s1gts2==false){

			i_save=i;
		
			s1gts2=true;
		}
		if((s1_current<=s2_current) && s1gts2==true)
		{
			//had to use buildvalue here to get rid of mem leak
			PyObject * current_pair;
			current_pair=PyList_New(2);
			PyList_SetItem(current_pair, 0,Py_BuildValue("O", PyList_GetItem(index,i_save))); 
			PyList_SetItem(current_pair, 1, Py_BuildValue("O", PyList_GetItem(index,i))); 

			PyList_Append(pairs,current_pair);
		
			s1gts2=false;
		}
	}
	
	if (s1gts2==true){

		PyObject * current_pair;
		current_pair=PyList_New(2);
		PyList_SetItem(current_pair, 0,Py_BuildValue("O",  PyList_GetItem(index,i_save))); 
		PyList_SetItem(current_pair, 1,Py_BuildValue("O",  PyList_GetItem(index,size-1))); 
		PyList_Append(pairs,current_pair);
		
	}

	return pairs;
error:
	return 0;
}


PyObject* threshold_trading(PyObject* self, PyObject* args)
{
	PyObject * price_series;

	PyObject * indicator_series;

	double sell_threshold;
	double buy_threshold;

	double initial_investment;
	double commission;

	PyObject* py_expectArgs;
	bool do_return_list;

	PyObject* py_expectArgs2;
	bool simulate;

	PyObject * return_list;

	if (!PyArg_ParseTuple(args, "O!O!ddddO!O!", &PyList_Type, &price_series, &PyList_Type, &indicator_series, &sell_threshold,&buy_threshold,&initial_investment,&commission, &PyBool_Type,&py_expectArgs, &PyBool_Type, &py_expectArgs2))
	{
		goto error;
	}

	int size;

	size = PyList_GET_SIZE(indicator_series);

	do_return_list= PyObject_IsTrue(py_expectArgs);

	simulate = PyObject_IsTrue(py_expectArgs2);

	return_list= PyList_New(size);

	bool baught;

	double baught_price;
	double investment = initial_investment;
	
	if (PyFloat_AsDouble(PyList_GetItem(indicator_series, 0))<= buy_threshold)
	{
		baught = true;
		baught_price = PyFloat_AsDouble(PyList_GetItem(price_series, 0));
		if (do_return_list)
		{
			PyList_SetItem(return_list, 0, PyFloat_FromDouble(investment/initial_investment));
		}
	}
	else
	{
		baught = false;
		baught_price = 0;
		if (do_return_list)
		{
			PyList_SetItem(return_list, 0, PyFloat_FromDouble(investment/initial_investment));
		}
	}



	for (int i = 1; i < size; i++)
	{
		double indicator_series_current;
		
		indicator_series_current = PyFloat_AsDouble(PyList_GetItem(indicator_series, i));
		


		if ((indicator_series_current >= sell_threshold) && baught == true)

		{
			//sell
			double price_series_current = PyFloat_AsDouble(PyList_GetItem(price_series, i));
			investment =  investment * price_series_current / baught_price - commission;
			baught = false;

		}

		if ((indicator_series_current <= buy_threshold) && baught == false)

		{
			//buy
			baught_price = PyFloat_AsDouble(PyList_GetItem(price_series, i));

			investment = investment - commission;
			baught = true;

		}
		if (simulate) { if (investment < 0) { investment = 0; } }
		if (do_return_list)
		{
			PyList_SetItem(return_list, i, PyFloat_FromDouble(investment/initial_investment));
		}
	}

	//PyObject * return_value=Py_flo
	if (do_return_list)
	{
		return return_list;
	}
	else
	{	
		return Py_BuildValue("d", investment / initial_investment);
	}
	
error:
	return 0;
}

PyObject* threshold_trading_binary(PyObject* self, PyObject* args)
{
	PyObject * price_series;

	PyObject * indicator_series;



	double initial_investment;
	double commission;

	PyObject* py_expectArgs;
	bool do_return_list;

	PyObject* py_expectArgs2;
	bool simulate;

	PyObject * return_list;

	if (!PyArg_ParseTuple(args, "O!O!ddO!O!", &PyList_Type, &price_series, &PyList_Type, &indicator_series, &initial_investment, &commission, &PyBool_Type, &py_expectArgs, &PyBool_Type, &py_expectArgs2))
	{
		goto error;
	}

	int size;

	size = PyList_GET_SIZE(indicator_series);

	do_return_list = PyObject_IsTrue(py_expectArgs);

	simulate = PyObject_IsTrue(py_expectArgs2);

	return_list = PyList_New(size);

	bool baught;

	double baught_price;
	double investment = initial_investment;

	if (PyFloat_AsDouble(PyList_GetItem(indicator_series, 0)) >0)
	{
		baught = true;
		baught_price = PyFloat_AsDouble(PyList_GetItem(price_series, 0));
		if (do_return_list)
		{
			PyList_SetItem(return_list, 0, PyFloat_FromDouble(investment / initial_investment));
		}
	}
	else
	{
		baught = false;
		baught_price = 0;
		if (do_return_list)
		{
			PyList_SetItem(return_list, 0, PyFloat_FromDouble(investment / initial_investment));
		}
	}



	for (int i = 1; i < size; i++)
	{
		double indicator_series_current;

		indicator_series_current = PyFloat_AsDouble(PyList_GetItem(indicator_series, i));



		if ((indicator_series_current <0) && baught == true)

		{
			//sell
			double price_series_current = PyFloat_AsDouble(PyList_GetItem(price_series, i));
			investment = investment * price_series_current / baught_price - commission;
			baught = false;

		}

		if ((indicator_series_current >0) && baught == false)

		{
			//buy
			baught_price = PyFloat_AsDouble(PyList_GetItem(price_series, i));

			investment = investment - commission;
			baught = true;

		}
		if (simulate) { if (investment < 0) { investment = 0; } }
		if (do_return_list)
		{
			PyList_SetItem(return_list, i, PyFloat_FromDouble(investment / initial_investment));
		}
	}

	//PyObject * return_value=Py_flo
	if (do_return_list)
	{
		return return_list;
	}
	else
	{
		return Py_BuildValue("d", investment / initial_investment);
	}

error:
	return 0;
}


PyObject* linear_regression(PyObject* self, PyObject* args)
{
	PyObject * values;
	PyObject * index;

	if (!PyArg_ParseTuple(args,"O!O!",&PyList_Type,&values,&PyList_Type,&index))
	{
		goto error;
	}

    int size;
	
	size=PyList_Size(values);

	PyObject *return_data = PyDict_New(); 

	double m;
	double b;
	double R;

	double x_squared_summation=0;
	double y_squared_summation=0;
	double xy_summation=0;
	double x_summation=0;
	double y_summation=0;

	for (int i=0;i<size;i++)
	{
		double current_index=PyFloat_AsDouble(PyList_GetItem(index,i));
		double current_value=PyFloat_AsDouble(PyList_GetItem(values,i));
		x_squared_summation = x_squared_summation + current_index*current_index;
		y_squared_summation = y_squared_summation + current_value*current_value;
		xy_summation = xy_summation + current_value*current_index;
		x_summation = x_summation + current_index;
		y_summation = y_summation + current_value;
	}


	double x_avg=x_summation /size;
	double y_avg=y_summation /size;


	m=(xy_summation-y_avg*x_summation-x_avg*y_summation+size*x_avg *y_avg )/(x_squared_summation-2*x_avg*x_summation+size*x_avg*x_avg);

	b=y_avg -m*x_avg;

	R=(xy_summation /size-x_avg *y_avg )/sqrt(((x_squared_summation /size)-x_avg *x_avg )*((y_squared_summation /size)-y_avg *y_avg) );
	R=R*R;

	return_data=Py_BuildValue ("{s:O,s:O,s:O}","m",Py_BuildValue("O",PyFloat_FromDouble(m)),"b",Py_BuildValue("O",PyFloat_FromDouble(b)),"R Squared",Py_BuildValue("O",PyFloat_FromDouble(R)));

	return return_data;
error:
	return 0;
}

PyObject* exponential_regression(PyObject* self, PyObject* args)
{
	PyObject * values;
	PyObject * index;

	if (!PyArg_ParseTuple(args,"O!O!",&PyList_Type,&values,&PyList_Type,&index))
	{
		goto error;
	}

    int size;
	
	size=PyList_Size(values);

	PyObject *return_data = PyList_New(3);

	double R;

	double x_squared_summation=0;
	double y_squared_summation=0;//add this to 
	double xy_summation=0;
	double x_summation=0;
	double y_summation=0;

	
	for (int i=0;i<size;i++)
	{

		double current_index=PyFloat_AsDouble(PyList_GetItem(index,i));

		double current_value=PyFloat_AsDouble(PyList_GetItem(values,i));

		x_squared_summation = x_squared_summation + current_index*current_index;

		y_squared_summation = y_squared_summation + log(current_value*current_value);

		xy_summation = xy_summation + log(current_value)*current_index;

		x_summation = x_summation + current_index;

		y_summation = y_summation + log(current_value);

	}


	double x_avg=x_summation /size;//mean x
	double y_avg=y_summation /size;

	double Sxx=x_squared_summation /size-x_avg *x_avg ;
	double Syy=y_squared_summation /size-(y_summation /size)*(y_summation /size);
	double Sxy=xy_summation/size-x_avg *y_avg ;


	double B=exp(Sxy /Sxx);
	double A=exp(y_avg -x_avg *log(B));


	R=Sxy /(sqrt(Sxx)*sqrt(Syy));
	R=R*R;


		//
	return_data=Py_BuildValue ("{s:O,s:O,s:O}","B (y=AB^x)",Py_BuildValue("O",PyFloat_FromDouble(B)),"A (y=AB^x)",Py_BuildValue("O",PyFloat_FromDouble(A)),"R Squared",Py_BuildValue("O",PyFloat_FromDouble(R)));

	return return_data;
error:
	return 0;
}

PyObject* polynomial_regression(PyObject* self, PyObject* args)
{
	//doesnt work

	PyObject * values;
	PyObject * index;
	PyObject * order;

	if (!PyArg_ParseTuple(args,"O!O!i",&PyList_Type,&values,&PyList_Type,&index,&order))
	{
		goto error;
	}

    int size;
	
	size=PyList_Size(values);

	PyObject *coefficients = PyList_New(PyFloat_AsDouble(order)+1);
	PyObject *return_data=PyDict_New ();

	double R;

	double x_squared_summation=0;
	double y_squared_summation=0;
	double xy_summation=0;
	double x_summation=0;
	double y_summation=0;
	
	
	for (int i=0;i<size;i++)
	{
		double current_index=PyFloat_AsDouble(PyList_GetItem(index,i));
		double current_value=PyFloat_AsDouble(PyList_GetItem(values,i));

		x_squared_summation = x_squared_summation + current_index*current_index;

		y_squared_summation = y_squared_summation + log(current_value*current_value);

		xy_summation = xy_summation + log(current_value)*current_index;
		x_summation = x_summation + current_index;

		y_summation = y_summation + log(current_value);

	}


	double x_avg=x_summation /size;//mean x
	double y_avg=y_summation /size;

	double Sxx=x_squared_summation /size-x_avg *x_avg ;
	double Syy=y_squared_summation /size-(y_summation /size)*(y_summation /size);
	double Sxy=xy_summation/size-x_avg *y_avg ;


	double B=exp(Sxy /Sxx);
	double A=exp(y_avg -x_avg *log(B));


	R=Sxy /(sqrt (Sxx)*sqrt(Syy));
	R=R*R;


		
	return_data=Py_BuildValue ("{s:O,s:O}","C0+C1x+C2x^2...",Py_BuildValue("O",coefficients),"R_squared",Py_BuildValue("O",PyFloat_FromDouble(R)));

	return return_data;
error:
	return 0;
}

PyObject* correlation_coefficient(PyObject* self, PyObject* args)
{
	PyObject * values;//values
	PyObject * index;//index that includes the dates

	if (!PyArg_ParseTuple(args,"O!O!",&PyList_Type,&values,&PyList_Type,&index))
	{
		goto error;
	}

    int size;//size of incoming list 
	
	size=PyList_Size(values);

	double m;
	double b;
	double R;

	double x_squared_summation=0;
	double y_squared_summation=0;
	double xy_summation=0;
	double x_summation=0;
	double y_summation=0;

	for (int i=0;i<size;i++)
	{
		double current_index=PyFloat_AsDouble(PyList_GetItem(index,i));
		double current_value=PyFloat_AsDouble(PyList_GetItem(values,i));
		x_squared_summation = x_squared_summation + current_index*current_index;
		y_squared_summation = y_squared_summation + current_value*current_value;
		xy_summation = xy_summation + current_value*current_index;
		x_summation = x_summation + current_index;
		y_summation = y_summation + current_value;
	}

	double x_avg=x_summation /size;
	double y_avg=y_summation /size;

	R=(xy_summation /size-x_avg *y_avg )/sqrt(((x_squared_summation /size)-x_avg *x_avg )*((y_squared_summation /size)-y_avg *y_avg) );
	//R=R*R;

	return Py_BuildValue("O",PyFloat_FromDouble(R));
error:
	return 0;
}

PyObject* list_dict_search(PyObject* self, PyObject* args)
{
	
	PyObject * list;

	long wo;

	if (!PyArg_ParseTuple(args,"O!l",&PyList_Type,&list,&wo))
	{
		goto error;
	}
	//return Py_BuildValue("i", wo);
    int size;
	
	size=PyList_Size(list);

	PyObject *return_data = PyList_New(0);
	
	bool found=false;
	
	for (int i=0;i<size;i++)
	{
		PyObject * dict=PyList_GetItem (list,i);
		
		if ( wo==PyLong_AsLong( PyDict_GetItemString(dict,"OrderNo"))){
			found=true;
			PyList_Append(return_data,Py_BuildValue("i", i));	
		}
		else{
			if (found==true){
				return return_data;
			}
		}
	}

	return return_data;
error:
	return 0;
}

PyObject* WOsearch(PyObject* self, PyObject* args)
{
	
	PyObject * list;

	double wo;

	if (!PyArg_ParseTuple(args,"O!d",&PyList_Type,&list,&wo))
	{
		goto error;
	}
	//return Py_BuildValue("i", wo);
    int size;
	
	size=PyList_Size(list);

	PyObject *return_data = PyList_New(0);
	
	
	for (int i=0;i<size;i++)
	{
		PyObject * current_line=PyList_GetItem (list,i);
		
		if ( wo==PyFloat_AsDouble(current_line)){

			PyList_Append(return_data,Py_BuildValue("i", i));	
		
		}
	}

	return return_data;
error:
	return 0;
}

PyObject* get_ema(PyObject* self, PyObject* args)
{
	PyObject * list;

	long period;

	if (!PyArg_ParseTuple(args,"O!l",&PyList_Type,&list,&period))
	{
		goto error;
	}

    int size;
	
	size=PyList_Size(list);

	PyObject *return_data = PyList_New(size);
	
	double multiplier=2./(period+1.);

	for (int i=0;i<period;i++)
	{
		PyList_SET_ITEM(return_data, i, PyFloat_FromDouble(0.)); 
	}

	for (int i=period;i<size;i++)
	{
		double ema=PyFloat_AsDouble (PyList_GetItem(return_data,i-1));

		PyObject *num =PyFloat_FromDouble((PyFloat_AsDouble( PyList_GetItem(list,i))-ema)*multiplier+ema);
	
		PyList_SET_ITEM(return_data, i, num); 
	}

	return return_data;
error:
	return 0;
}


PyObject* test(PyObject* self, PyObject* args)
{
	PyObject * TP;

	

	PyObject * RMF;

	int period=14;

	if (!PyArg_ParseTuple(args,"O!O!",&PyList_Type,&TP,&PyList_Type,&RMF))
	{
		goto error;
	}

    int size;
	size=PyList_GET_SIZE(TP);

	PyObject *lst = PyList_New(size);
	//allocate the array
	for (int i=0;i<period;i++){
		
		PyObject *num = PyFloat_FromDouble(0);
		PyList_SET_ITEM(lst, i, num); 
	}

	for (int i=period;i<size;i++){
		double pos=0;
		double neg=0;	
		for (int k = i-period+1;k<i;k++){		
			double misc=PyFloat_AsDouble(PyList_GetItem(TP,k));
			double misc2=PyFloat_AsDouble(PyList_GetItem(TP,k-1));
		
			if (misc>misc2){
				pos=pos+ PyFloat_AsDouble( PyList_GetItem(RMF,k));
			}
			else
			{
				neg=neg+ PyFloat_AsDouble( PyList_GetItem(RMF,k));
			}	
		}
		PyObject *num = PyFloat_FromDouble(100-100/(1+pos/neg));
		PyList_SET_ITEM(lst, i, num); 		
	}

	return lst;
error:
	return 0;
}


//

PyDoc_STRVAR(william_r__doc__,"x,y,max_iterations");
PyDoc_STRVAR(tech_analysis_c__doc__,"Implements C funcions for Securities Wizard");

static PyMethodDef SpamMethods[] = {
	{"william_r",(PyCFunction)william_r,METH_VARARGS,william_r__doc__},
	{"money_flow_index",(PyCFunction)money_flow_index,METH_VARARGS,0},	
	{"ultimate_ocillator",(PyCFunction)ultimate_ocillator,METH_VARARGS,0},
	{"find_spans_one_above_other",(PyCFunction) find_spans_one_above_other,METH_VARARGS,0},
	{"linear_regression",(PyCFunction) linear_regression,METH_VARARGS,0},
	{"exponential_regression",(PyCFunction) exponential_regression,METH_VARARGS,0},
	{"polynomial_regression",(PyCFunction) polynomial_regression,METH_VARARGS,0},
	{"correlation_coefficient",(PyCFunction) correlation_coefficient,METH_VARARGS,0},
	{"threshold_trading",(PyCFunction)threshold_trading,METH_VARARGS,0},
	{"threshold_trading_binary",(PyCFunction)threshold_trading_binary,METH_VARARGS,0},
	{"list_dict_search",(PyCFunction) list_dict_search,METH_VARARGS,0},
	{"WOsearch",(PyCFunction) WOsearch,METH_VARARGS,0},
	{"get_ema",(PyCFunction) get_ema,METH_VARARGS,0},
	{"test",(PyCFunction) test,METH_VARARGS,0},
    {0, 0, 0, 0}        /* Sentinel */
};

static struct PyModuleDef tech_analysis_cmodule = {
   PyModuleDef_HEAD_INIT,
   "tech_analysis_c",   /* name of module */
   tech_analysis_c__doc__, /* module documentation, may be NULL */
   -1,       /* size of per-interpreter state of the module,
                or -1 if the module keeps state in global variables. */
   SpamMethods
};

PyMODINIT_FUNC PyInit_tech_analysis_c(void)
{
    PyObject *m;

    m = PyModule_Create(&tech_analysis_cmodule);
    if (m == NULL)
        return NULL;

    
    return m;
}