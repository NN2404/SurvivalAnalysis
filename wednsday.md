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
 
 
  
     