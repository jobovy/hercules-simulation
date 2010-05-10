PYTHON=python
BARDIR=../bar

bar: $(BARDIR)/rphi2d.ps $(BARDIR)/rphi1d.ps $(BARDIR)/detecta.ps \
	$(BARDIR)/detectb.ps

$BARDIR/rphi2d.ps: bar_figures.py calc_veldist_2d.py
	$(PYTHON) bar_figures $@

$BARDIR/rphi1d.ps: bar_figures.py calc_veldist_1d.py
	$(PYTHON) bar_figures $@ vlos

$(BARDIR)/detecta.ps: bar_detectability.py
	$(PYTHON) bar_detectability.py $@

$(BARDIR)/detectb.ps: bar_detectability.py
	$(PYTHON) bar_detectability.py $@ --skipCenter=0.15

