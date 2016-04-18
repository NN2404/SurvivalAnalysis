Survival analysis is the phrase used to describe the analysis of data in the form of times from a well-defined time origin until the 
occurrence of some particular event or end point. The time origin refers to the time in which an individual enters an experimental 
study EG.A clinical trial comparing 2 treatments. This time maybe measured in days, months or years depending on the experiment. 
The event might refer to the death of a patient, then the recorded data is literally survival times, 
but the event isn’t limited to the death of a patient. In some cases it might refer to the recurrence of symptoms or in
the case of motor vehicles it might refer to the length of time a particular component will last. 

<pre><code>
 
loglikePHexp <- function(theta,data){

    ti <- data[,1]
    di <- data[,2]
    xi <- as.matrix(data[,-(1:2)])


    lam <- exp(theta[1])
    betas <- theta[-1]
    phi <- exp(xi%*%betas)

    loglike <- sum(di*(log(phi) + log(lam)) - phi*lam*ti)
    -loglike
}

 

</code></pre>




Survival data is generally not normally distributed. A histogram made up of survival times of a group in an experiment
will usually form a positively skewed histogram which will contain a long “tail” to the right of the data.This is one 
reason why survival data is unsuitable to the regular statistical procedures used in survival analysis. The main reason
is that survival times are nearly always censored.

The survival time of a particular individual is said to be censored the end-point of interest has not been observed for
that individual. This may occur for a number of reasons. If an individual in the experiment decides to leave the experiment
before the event of interest has occurred for that individual, their data is said to be censored.  
 
Censoring in Survival analysis
The survival time of a particular individual is said to be censored the end-point of interest has not been observed for that individual. This may occur for a number of reasons. If an individual in the experiment decides to leave the experiment before the event of interest has occurred for that individual, their data is said to be censored. Also if the experiment trial finishes before the event occurs for a particular individual, their data is also said to be censored. An actual survival time is also said to be censored if an individual in the experiment dies of a cause unrelated to the treatment they were on eg. traffic accident.
There are 3 types of censoring; Right Censoring, Left censoring, and Interval Censoring.
Right Censoring:
In all survival analysis cases the individual arrives into the study at time t0 and they die at time t0 + t. In uncensored data t is known whereas in censored data t is unknown, either because a person was lost to follow-up or they were alive a when the experiment finished. If the latest a person was known to be alive at t0 + c then c is that person’s censored survival time.  If you were to illustrate this person’s data on a time line the event of interest for this persons occurs to the right of this persons last known survival time and so is called right-censoring. In short it means the event occurred after a certain point but it is unknown by how much it exceeds this point.
Left Censoring
To explain this form of censoring, consider an experimental study where the end point of interest is the recurrence of symptoms of a heart condition after someone undergoes heart surgery. One month after the operation, each of the individuals in the study is examined to see if any of the symptoms has reoccurred. During these tests some of the pre-operation symptoms have reoccurred for one of the individuals in the study. For this individual the reoccurrence time is less than one month. It is impossible to isolate the exact time that the reoccurrence took place. Again if you were to illustrate this person’s survival data on a time line the actual time of reoccurrence is to the left  of the person’s check-up and so is called left censoring.  In other words the event took place before a known time but it is unknown by how much.  Left censoring is far less common the right censoring.
Interval Censoring
Another type of censoring is interval censoring. This is where an event is known to have occurred for an individual during a known period of time, but it is not known the exact time this event occurred. Again consider the example for left censoring. During the first check-up after one month a certain individual is known to be free of symptoms one but during their second check-up after two months that individual starts to experience recurring symptoms. It is impossible to determine the examine the exact time in which symptoms started reoccurring, all we know is that it happened between the individuals first and second check-up. This individual’s statistical data would therefore have to be censored and so this form of censoring is known as interval censoring.  Interval censoring incorporates both right censoring. In essence it is censoring between two points in time (say t1 and t2). If t1 is 0 then it is simply left censoring.
By far the most common form of censoring is right censoring and that will be the main form of censoring that I will be covering during this project.

  The survival function, also known as a survivor function or reliability function, is a property of any random variable that maps a set of events, usually associated with mortality or failure of some system, onto time. It captures the probability that the system will survive beyond a specified time.

The term reliability function is common in engineering while the term survival function is used in a broader range of applications, including human mortality. Another name for the survival function is the complementary cumulative distribution function.
   
Let T be a continuous random variable with cumulative distribution function F(t) on the interval [0,∞). Its survival function or reliability function is:
<pre><code>
   S(t) = int_{t}^{\infinite} f(u) du
  
</code></pre>


Every survival function S(t) is monotonically decreasing. The time, t = 0, represents some origin, typically the beginning of a study or the start of operation of some system. S(0) is commonly unity but can be less to represent the probability that the system fails immediately upon operation. Since the CDF is a right-continuous function, it can be said that the survival function S(t) = 1-F(t) is also right-continuous.
