<TeXmacs|1.0.7.14>

<style|generic>

<\body>
  Initially the model parameters are scaled by choosing
  <math|\<alpha\><rsub|j>> such that <math|\<theta\><rsub|j>/\<alpha\><rsub|j>=<wide|\<theta\>|^><rsub|j>\<in\><around*|[|0,1|]>>
  and the vector of scaled parameters defined as
  <math|<with|math-font-series|bold|<wide|\<theta\>|^>>=<around|(|<wide|\<theta\>|^><rsub|1>,\<ldots\>,<wide|\<theta\>|^><rsub|p>|)>\<in\><R><rsup|p>>.
  Considering the specific response <math|x<rsub|i>> we then have,

  <\equation>
    <label|eq:asf-scale>\<alpha\><rsub|j>*<frac|\<partial\>*x<rsub|i>|\<partial\>*\<theta\><rsub|j>>=<frac|\<partial\>*x<rsub|i>|\<partial\>*<wide|\<theta\>|^><rsub|j>>,
  </equation>

  which is an important result when comparing to other sensitivity measures.
  Taking a Taylor expansion of <math|x<rsub|i><around|(|t,<with|math-font-series|bold|<wide|\<theta\>|^>>|)>>
  about the point <math|<wide|\<theta\>|^><rsub|j>\<pm\>\<Delta\>*<wide|\<theta\>|^><rsub|j>>,
  where we have abused notation slightly defining
  <math|x<rsub|i><around|(|t,<with|math-font-series|bold|\<theta\>>|)>=x<rsub|i><around|(|t,<with|math-font-series|bold|<wide|\<theta\>|^>>|)>>,
  one obtains,

  <\equation>
    <label|eq:asf-expand><frac|\<partial\>*x<rsub|i><around|(|t,<with|math-font-series|bold|<wide|\<theta\>|^>>|)>|\<partial\>*<wide|\<theta\>|^><rsub|j>>\<approx\><frac|x<rsub|i><around|(|t,<with|math-font-series|bold|\<Delta\>*<wide|\<theta\>|^>>|)>-x<rsub|i><around|(|t,<with|math-font-series|bold|<wide|\<theta\>|^>>|)>|\<pm\>\<Delta\>*<wide|\<theta\>|^><rsub|j>>=E<rsub|j><around|(|<math-bf|x>|)>,
  </equation>

  where <math|<with|math-font-series|bold|\<Delta\>*<wide|\<theta\>|^>>=<around|(|<wide|\<theta\>|^><rsub|1>,<wide|\<theta\>|^><rsub|2>,\<ldots\>,<wide|\<theta\>|^><rsub|j>\<pm\>\<Delta\>*<wide|\<theta\>|^><rsub|j>,\<ldots\>,<wide|\<theta\>|^><rsub|p>|)>>.
  Essentially the parameters <math|\<theta\><rsub|j>> are varied
  one-at-a-time (OAT) to obtain an approximation of the local gradient of the
  function <math|<math-bf|x>>. The scaled quantity,
  <math|E<rsub|j><around|(|<math-bf|x>|)>>, related to the absolute
  sensitivity function, is defined as the elementary effect of the <math|j>th
  factor at a point <math|<with|math-font-series|bold|<wide|\<theta\>|^>>>.
</body>

<\initial>
  <\collection>
    <associate|language|american>
    <associate|preamble|false>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|eq:asf-expand|<tuple|2|?>>
    <associate|eq:asf-scale|<tuple|1|?>>
  </collection>
</references>