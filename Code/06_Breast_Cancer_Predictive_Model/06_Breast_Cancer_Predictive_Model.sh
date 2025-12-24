#1 100-fold repeated CV model predictions + running
for i in {1..10}; #10 folds
do for l in {1..100}; #100 repeats
do sed  s/foldno/${i}/g breast.cancer.discovery.cv.model.R  | sed  s/seedno/${l}/g   > breast.cancer.discovery.cv.${i}.${l}.R;
./breast.cancer.discovery.cv.${i}.${l}.R;
done; done

#2 100-fold CV performance discovery
./breast.cancer.discovery.cv.performance.R

#test set model + performance
./breast.cancer.testset.R

#external data model + performance
./breast.cancer.testset.cv.model.R
