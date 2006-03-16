/* shouldReadFromFiles(my_sem_t)
 *  caso o count for 0, devolver true.
 *  cc, devolver false
 *  apos o teste e antes de retornar, incrementar o count.
 */
int shouldRestoreFromFiles(my_sem_t *sem)
{

	pthread_mutex_lock(&sem->lock);
	if(sem->count++ == 0) {
		pthread_mutex_unlock(&sem->lock);
		return 1;
	} else {
		pthread_mutex_unlock(&sem->lock);
		return 0;
	}
}

/** shouldWriteToFiles(my_sem_t, int max)
 * decrementar o count.
 * caso o count for 0, devolver true.
 * cc, devolver false.
 */
int shouldSaveToFiles(my_sem_t *sem)
{

	pthread_mutex_lock(&sem->lock);
	sem->count--;
	if(sem->count == 0) {
		pthread_mutex_unlock(&sem->lock);
		return 1;
	} else {
		pthread_mutex_unlock(&sem->lock);
		return 0;
	}
}

char* DumpCorrData(struct aiGridData *aiData)
{
	//TODO: mudar isto para nao usar tantos arrays (2x1024)
	char buf[1024];
	char ret[1024];
	int i;

	strcpy(ret, "");
	for(i=0; i < aiData->layersNum; i++) {	
		// printf("%d of %d: %Lf\n", i, this->layers_num, this->cr_med[i]);
		strcpy(buf, "");
		sprintf(buf, "\t\t\tLayer %2d: %.5f\n", i+1, aiData->layersCorr[i]);
		strcat(ret, buf);
	}
	strcat(ret, "\t\t\t-------------------\n");
	sprintf(buf, "\t\t\t  Global: %.5f\n", aiData->globalCorr);
	strcat(ret, buf);

	return strdup(ret);
}

