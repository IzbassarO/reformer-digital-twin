# SCOPE & LITERATURE BRIEF: Physics-informed digital twin of an industrial SMR (process + TWT + creep-life under uncertainty)

## Резюме (executive summary)
Ниша подтверждена как **ОТКРЫТАЯ**. Ни одна опубликованная работа не объединяет в одной рамке (i) 1D process-модель/ML-суррогат SMR, (ii) прогноз tube-wall temperature (TWT), (iii) оценку creep-life / remaining-life труб (HK40/HP-Nb) и (iv) количественную оценку неопределённости (UQ) на исключительно открытых данных. Литература распадается на три «силоса»: (а) process/furnace↔TWT (Plehiers & Froment 1989; Zamaniyan 2008; Latham 2011; Kumar/Baldea/Edgar 2015–2017; UCLA-группа Christofides 2016–2021), (б) operation↔tube life (Applied Sciences 2021 — ближайший аналог, детерминированный CFD; ERA «REFORM» — вероятностный, но проприетарный), (в) creep↔uncertainty без process-модели (Ray 2003; IJPVP 2023 neuro-fuzzy θ; IJPVP 2025 HP40Nb ML). Группа Christofides решает температурную равномерность и максимизацию H₂, но НЕ ресурс труб. Рекомендуемая конфигурация: 1D псевдогомогенная модель (Xu–Froment) + упрощённая furnace-модель (heat-flux/зонная по Latham) + LMP/Robinson creep-модуль + суррогат (gradient boosting/GP) + UQ (Sobol/SALib, conformal/MAPIE, NSGA-II/pymoo). Открытые creep-данные: NIMS CDS 16B (HK40) и 38A (HP40). Целевой журнал IJHE: 8000 слов, ≤12 рисунков, 3–5 highlights, structured abstract 200 слов, обязательная декларация об ИИ.

## 1. Проверка ниши (gap verification)

**Вывод: ниша ПОДТВЕРЖДЕНА КАК ОТКРЫТАЯ.** Существуют попарные связки, но нет ни одной работы со всеми четырьмя элементами на открытых данных. Отраслевая значимость темы высока: по данным Sweeney et al. 2023 (Ind. Eng. Chem. Res., 10.1021/acs.iecr.3c02045), «Steam methane reforming (SMR) currently supplies 76% of the world's hydrogen (H2) demand, totaling ∼70 million tonnes per year»; при этом фактический ресурс труб сильно варьируется — по Ray et al. 2003 (Eng. Fail. Anal. 10:351–362) «The actual service life varies from 30,000 to 180,000 [h] depending on the quality of materials and the service conditions», а преждевременные разрушения в первые 3–8 лет распространены. Это прямо мотивирует tube-life-aware эксплуатацию.

| Работа | DOI | Что связано | Данные | Отличие нашей работы |
|---|---|---|---|---|
| Plehiers & Froment 1989, Chem. Eng. Technol. 12:20–26 | 10.1002/ceat.270120105 | Зонный теплообмен топки + 1D гетерогенный реактор → температурное поле, эффлюент | Промышленная установка | Нет creep/life, нет UQ, нет ML |
| Zamaniyan, Ebrahimi, Soltan Mohammadzadeh 2008, Chem. Eng. Process. 47:946–956 | 10.1016/j.cep.2007.03.005 | 1D гетер. реактор + 3D зонная топка → внешняя TWT (max ~1/3 длины сверху) | Промышленный top-fired | Нет creep-life, нет UQ |
| Latham, McAuley, Peppley, Raybold 2011, Fuel Process. Technol. 92:1574–1586 | 10.1016/j.fuproc.2011.04.001 | Зонный метод Hottel; профили TWT (внешн./внутр.), furnace/process gas | Промышленные данные, калибровка | Мониторинг TWT, но НЕ creep/remaining-life |
| Kumar, Baldea, Edgar, Ezekoye 2015, IECR 54:4360–4370 | 10.1021/ie504087z | Smart-manufacturing: сенсоры+soft sensor+ROM/CFD, балансировка TWT | Промышленный SMR | Цель — энергоэффективность/равномерность, не ресурс |
| Kumar, Baldea, Edgar 2017, Comput. Chem. Eng. 105:224–236 (physics-based) | 10.1016/j.compchemeng.2017.01.002 | Physics-based furnace, неоднородное T-поле | — | Нет creep-life |
| Lao, Aguirre, Tran, Wu, Durand, Christofides 2016, Chem. Eng. Sci. 148:78–92 | 10.1016/j.ces.2016.03.038 | CFD одиночной трубы + управление; ограничение по max OTWT | CFD | Только температура/управление, НЕ ресурс труб |
| Tran, Aguirre, Durand, Crose, Christofides 2017, Comput. Chem. Eng. 104:185–200 | 10.1016/j.compchemeng.2017.04.026 | CFD/data-based балансировка топки | CFD | Равномерность TWT, не creep |
| Tran, Aguirre, Durand, Crose, Christofides 2017, Chem. Eng. Sci. 171:576–598 | 10.1016/j.ces.2017.06.001 | 3D CFD промышленной топки | CFD | Тяжёлый CFD; не creep; не для ноутбука |
| «Burner operation & catalyst tube lifetime», Applied Sciences 2021, 11:231 | 10.3390/app11010231 | operation→TWT→tube life (Larson-Miller) — **БЛИЖАЙШИЙ аналог** | Детерминированный CFD | Наш: 1D+ML+UQ, без больших CFD-кампаний; добавляем UQ и Pareto |
| Latham et al. 2017, Chem. Eng. Res. Des. (online TSK) S0263876217304173 | 10.1016/j.cherd.2017.09.023 | Онлайн-прогноз tube skin T; упоминает «digital twin» топки | Промышленный top-fired | Нет интегрированного creep-damage под неопределённостью |
| Ray, Sinha, Tiwari, Swaminathan et al. 2003, Eng. Fail. Anal. 10:351–362 | 10.1016/S1350-6307(02)00029-8 | Анализ разрушения труб; ускоренный creep, LMP | Служебные трубы (HK-40 mod.) | Материаловедение без process-модели |
| Neuro-fuzzy θ-projection HP40, IJPVP 2023 | S0308016123000558 (10.1016/j.ijpvp.2023.104909) | θ-projection + neuro-fuzzy, разброс creep-параметров | Служебная труба HP40, 870 °C | Creep-модуль без связи с process/TWT |
| HP40Nb ML microstructure, IJPVP 2025, vol. 216 | S0308016125000961 | ML-информированная эволюция микроструктуры, 950 °C/1600 ч | Экспозиция HP40Nb | Микроструктура, не системный digital twin |
| ERA «REFORM» (Brear et al. 2001, IJPVP 78:985–994) | 10.1016/S0308-0161(01)00113-2 | Вероятностный creep+thermal cycling, Монте-Карло | Проприетарно, НЕ открыто | Наш — на открытых данных, с process/ML-связкой |

**UCLA/Christofides:** группа (Lao, Aguirre, Tran, Wu, Crose, Christofides, 2016–2019; плюс IECR 2020 10.1021/acs.iecr.0c00456) решает температурную равномерность (temperature balancing), балансировку топки и максимизацию производства H₂ под ограничением максимальной OTWT ≤ проектной. Ресурс труб / накопление ползучести / remaining-life они **НЕ моделируют** — только держат TWT ниже проектного лимита как жёсткое ограничение. Это ключевое отличие нашей работы.

## 2. Состояние по модулям + рекомендации

**(a) Реактор.** Варианты: 1D псевдогомогенная / 1D гетерогенная (с фактором эффективности) / 2D. Кинетика: Xu & Froment 1989 (Ni/MgAl₂O₄, LHHW) — де-факто стандарт; альтернативы Hou & Hughes 2001, Numaguchi & Kikuchi 1988. Факторы эффективности промышленных Ni-катализаторов малы (0.01–0.1) из-за диффузионных ограничений. Деактивация (спекание, углеотложение, сера) — через множитель активности a(t)∈(0,1]. **Рекомендация:** 1D псевдогомогенная модель с кинетикой Xu–Froment и постоянным/скалярным фактором эффективности + перепад давления по Эргуну. Достаточно для промышленной валидации, принимается рецензентами, считается на ноутбуке.

**(b) Топка/теплообмен.** Top-fired vs side-fired; методы: зонный (Hottel), Roesler flux, заданный профиль теплового потока, «упрощённые furnace-модели» (Plehiers & Froment 1989; Latham 2011; Kumar–Baldea–Edgar 2015–2017; Pantoleontos 2012; Zamaniyan 2008; Pedernera 2003; Rajesh 2000; Nandasana 2003). TWT-профили прогнозируют Latham 2011 (зонный Hottel, валидирован на промышленных данных, max TWT ~1/3 сверху) и Zamaniyan 2008 (3D-зонный). **Рекомендация:** заданный/параметризованный профиль теплового потока q(z) ИЛИ упрощённая 1D зонная модель типа Latham; это даёт правдоподобный профиль внешней TWT без больших CFD. CFD-кампании исключить (несовместимы с Colab/ноутбуком).

**(c) TWT.** Цепочка сопротивлений: газ топки → (радиация+конвекция) → внешняя стенка → теплопроводность стенки → внутренняя стенка → процессный газ. Рабочая температура стенки обычно 850–1095 °C, а локально может достигать 1150–1200 °C (US Patent RE50475: «the reformer tubes operate at a highly elevated temperature (up to 1150-1200° C.) such that they are susceptible to … "creep"»). IR-пирометрия / «tube skin temperature survey» дают «hot bands», flame impingement. Потеря активности катализатора → рост TWT (эндотерма не поглощает тепло). Основа проектирования — API 530. **Рекомендация:** аналитическая радиально-теплопроводная цепочка + q(z) из furnace-модуля; hot bands задавать как локальные возмущения профиля.

**(d) Creep-life труб.** Сплавы: HK40 (25Cr-20Ni-0.4C, SCH22-CF) и HP/HP-Nb (25Cr-35Ni-Nb, SCH24, «micro-alloyed»). Методы: Larson–Miller (LMP = T(C + log t_r), C≈15–20; API 530 / данные производителей), θ-projection, правило Robinson (life-fraction) для переменных T/σ, hoop stress по формуле среднего диаметра. Трубы проектируют по API STD 530 на средний ресурс до разрушения по ползучести 100 000 ч (≈11,4 лет) (Heat Exchanger World / ScienceDirect Topics: «Reformer tubes are normally designed according to the "remaining life assessment" technique, API-530, for an average lifetime before creep rupture of 100,000 hours»). Правило чувствительности: согласно US Patent 11215574, «the life time will typically be reduced to 50% if the operating temperature consistently is 10 to 15° C. higher» — то есть распространённое эвристическое «+20 °C → ×0.5» является менее консервативным, чем цитируемый диапазон 10–15 °C; в статье лучше приводить обе оценки со ссылками.

Инвентарь открытых creep-данных:

| Источник | DOI/PII | Сплав | Что в таблицах | Доступ |
|---|---|---|---|---|
| NIMS CDS 16B | 10.11503/nims.1020 | HK40 (SCH22-CF, 25Cr-20Ni-0.4C) | Четырнадцать плавок центробежнолитой трубы; stress-rupture до ~100 000 ч при 800–1000 °C (≈263 точек t_r по вторичным источникам) | Регистрация NIMS; условия использования |
| NIMS CDS 38A | 10.11503/nims.1042 | HP40 (SCH24, 25Cr-35Ni-0.4C) | Многоплавочные stress-rupture данные | NIMS, условия использования |
| Wilshire & Whittaker 2013, MSEA | S0921509313005716 | HK40, HP40 | Нормализованные кривые на базе NIMS; в осн. графики | Подписка |
| Ray et al. 2003, EFA 10:351–362 | 10.1016/S1350-6307(02)00029-8 | HK-40 mod. | Ускоренный creep, LMP, hot tensile | Подписка |
| Neuro-fuzzy θ HP40, IJPVP 2023 | S0308016123000558 | HP40 (870 °C, 50–68 МПа) | θ-параметры, разброс; частично графики | Подписка |
| HP40Nb ML, IJPVP 2025, v.216 | S0308016125000961 | HP40Nb (950 °C/1600 ч) | Доли фаз, нанотвёрдость; в осн. графики | Подписка |
| Ghatak & Robi, HP40Nb creep 950 °C, JMST 2019 | 10.1007/s12206-019-0922-9 | HP40Nb | Кривые ползучести, LMP | Подписка |
| Le May, da Silveira, Vianna 1996, IJPVP 66:233–241 | 10.1016/0308-0161(95)00098-4 | Reformer tubes | Критерии повреждения/ресурса | Подписка |
| Alvino et al. 2010, EFA 17:1526–1541 | 10.1016/j.engfailanal.2010.06.003 | HP tubes ~10 лет | Микроструктура/повреждение | Подписка |
| European J. Materials 2021 (short-term overheating HP-Nb) | 10.1080/26889277.2021.1994841 | HP-Nb | Creep-rupture при кратком перегреве | Open Access (CC BY) |

Даташиты производителей (Schmidt+Clemens Centralloy G4852 Micro/ET45 Micro; Kubota KHR35C/KHR45A; Manoir Manaurite XM/XTM) содержат stress-rupture кривые (обычно графики → нужна оцифровка). **Рекомендация:** базой сделать NIMS 16B/38A (табличные, многоплавочные — идеальны для heat-to-heat разброса), калибровать LMP; European J. Materials 2021 (OA) — для сценариев перегрева.

**(e) ML-суррогаты + UQ.** Дизайн: Latin hypercube (LHS), обычно 300–2000 точек. Модели: gradient boosting (XGBoost/LightGBM) и Gaussian processes для UQ; MLP как альтернатива. UQ: conformal prediction (MAPIE), quantile regression, deep ensembles. Метрики: R², RMSE, PICP, MPIW. Глобальная чувствительность: Sobol (SALib). Многоцелевая оптимизация: NSGA-II (pymoo). Методические шаблоны: arXiv:2507.07641 (Nabavi, Guo, Wang 2025 — ML-суррогат SMR + NSGA-II MOO); Jafarizadeh, Panjepour, Davazdah Emami 2024, IJHE 96:1262–1280 (10.1016/j.ijhydene.2024.11.352 — CFD→ML(DNN)→SHAP→оптимизация); Yerimah, Mehana, Singh, Sharan 2025, IECR 64:18298–18314 (10.1021/acs.iecr.5c01681 — «Uncertainty-Aware» оптимизация SMR+WGS+PSA).

**(f) ML для creep-life жаропрочных сплавов (2019–2026).** Методы: ANN, XGBoost/gradient boosting, GP, transfer learning; данные — преимущественно NIMS. Heat-to-heat разброс — через включение состава/микроструктуры как признаков и вероятностные модели. Physics-informed признаки: LMP, Wilshire. Примеры: ISIJ Int. 63(10) 2023 (ML на NIMS для ферритных сталей, GP/NN/GBDT); IJPVP 2019 (ANN для 9Cr-1Mo-V-Nb, S0308016119302844); IJPVP 2024 (probabilistic creep, superheater header, material uncertainty, 10.1016/j.ijpvp.2024.105211). **Рекомендация:** не переизобретать общий ML-creep; наш вклад — связка creep-модуля с process/TWT и распространение неопределённости на remaining-life.

**(g) Гибкая/динамическая эксплуатация H₂-установок.** SMR-заводы обычно безопасно наращивают нагрузку ~10% от проектной за интервал; риск «over firing» при быстрых скачках. Транзиенты → термическая усталость + creep-fatigue; релаксация термических напряжений даёт основную долю повреждения при переходных режимах (NACE CORROSION 1998: при неучёте транзиентов оценка ресурса нонконсервативна в ≥10 раз, а термические напряжения могут на порядок превышать hoop stress от давления). Сопряжение SMR с электролизом/ВИЭ повышает переменность. **Рекомендация:** использовать для обоснования сценариев RQ2.

## 3. Инвентарь валидационных кейсов

| Кейс | DOI/источник | Что сообщается | Таблицы/рисунки |
|---|---|---|---|
| Xu & Froment 1989 Part II | 10.1002/aic.690350110 | Промышленный реформер: профили конверсии/парц. давлений, факторы эффективности; типовые d_i≈0.1 м, L≈12 м, S/C~3 | Профили — рисунки; часть параметров — текст; **точные табличные значения не подтверждены (пейволл)** |
| Latham 2011 | 10.1016/j.fuproc.2011.04.001 | Профили внешней/внутренней TWT, furnace/process gas T; геометрия, состав; подгонка к промышленным данным | Смешанно; профили TWT валидированы |
| Latham 2008 (диссертация Queen's) | queensu.scholaris.ca | Полное описание геометрии топки, S/C 1–4, подогрев процесса ~565 °C | Диссертация — детально |
| Pantoleontos, Kikkinides, Georgiadis 2012, IJHE 37:16346–16358 | 10.1016/j.ijhydene.2012.02.125 | Гетер. динамич. модель; заданный квадратичный профиль T стенки; оптимизация H₂ | Профили — рисунки |
| Kumar, Baldea, Edgar, Ezekoye 2015, IECR 54:4360–4370 | 10.1021/ie504087z | 336 труб (7×48), 96 горелок (8×12); распределение max TWT | Рисунки + параметры |
| Zamaniyan 2008, Chem. Eng. Process. 47:946–956 | 10.1016/j.cep.2007.03.005 | 3D-зонная топка, внешняя TWT профиль | Профили — рисунки |
| Olivieri & Vegliò 2008, Fuel Process. Technol. | (DOI не извлечён) | Оптимизация распределения топлива в топке | — |
| Plehiers & Froment 1989 | 10.1002/ceat.270120105 | Температурное поле топки + эффлюент | Рисунки |
| Lao 2016 / Tran 2017 (UCLA CFD) | 10.1016/j.ces.2016.03.038 и др. | Детальные CFD-профили TWT | Рисунки; тяжёлый CFD |

**Топ-3 для валидации:**
1. **Latham 2011 + диссертация Latham 2008** — единственный полностью документированный кейс с валидированными профилями внешней TWT; недостаёт: точных heat-flux данных по зонам.
2. **Xu & Froment 1989 Part II** — канонический кейс реактора (конверсия/состав); недостаёт: подтверждённых табличных TWT/heat-flux (нужен полный PDF).
3. **Pantoleontos 2012 (IJHE)** — соответствует целевому журналу, есть заданный профиль T стенки и валидация против литературы; недостаёт: измеренной TWT (профиль задан, не измерен).

## 4. Верифицированная библиография (DOI)
- Xu & Froment 1989 I, AIChE J. 35(1):88–96 — 10.1002/aic.690350109 ✓
- Xu & Froment 1989 II, AIChE J. 35(1):97–103 — 10.1002/aic.690350110 ✓
- Plehiers & Froment 1989, Chem. Eng. Technol. 12:20–26 — 10.1002/ceat.270120105 ✓
- Latham et al. 2011, Fuel Process. Technol. 92:1574–1586 — 10.1016/j.fuproc.2011.04.001 ✓
- Pantoleontos et al. 2012, IJHE 37:16346–16358 — 10.1016/j.ijhydene.2012.02.125 ✓
- Kumar, Baldea, Edgar, Ezekoye 2015, IECR 54(16):4360–4370 — 10.1021/ie504087z ✓
- Kumar, Baldea, Edgar 2016, Control Eng. Practice 54:140–153 — 10.1016/j.conengprac.2016.05.010 (по вторичным источникам — сверить)
- Zamaniyan et al. 2008 — **ИСПРАВЛЕНО:** вероятно Chem. Eng. Process. 47(5):946–956, 10.1016/j.cep.2007.03.005 (НЕ J. Nat. Gas Chem.; исходная атрибуция ошибочна)
- Wilshire & Whittaker 2013, MSEA — S0921509313005716 ✓
- arXiv:2507.07641 (Nabavi, Guo, Wang 2025 — принята на IEEE QSW 2025) ✓
- IECR 2025 10.1021/acs.iecr.5c01681 (Yerimah, Mehana, Singh, Sharan; 64(37):18298–18314) ✓
- IJHE 2024 CFD→ML→SHAP, S0360319924050432 (Jafarizadeh et al., 96:1262–1280) — 10.1016/j.ijhydene.2024.11.352 ✓
- NIMS CDS 16B — 10.11503/nims.1020 ✓; NIMS CDS 38A — 10.11503/nims.1042 ✓
- Ray et al. 2003, EFA 10:351–362 — 10.1016/S1350-6307(02)00029-8 ✓
- European J. Materials 2021 — 10.1080/26889277.2021.1994841 ✓ (**НЕ Materials at High Temperatures**; журнал — European Journal of Materials; тема — creep rupture in HP-Nb tubes due to short-term overheating)
- IJPVP 2023 neuro-fuzzy θ HP40 (Guguloth & Roy) — S0308016123000558 ✓
- IJPVP 2025 HP40Nb ML (Khan, Shunmugasamy, … Laycock, Mansoor), v.216 — S0308016125000961 ✓

## 5. Требования журнала
**IJHE (2026):** Research Papers ≤8000 слов (включая таблицы, без abstract/references/SI), ≤12 рисунков; Review ≤12000 слов; Short Communication ≤3000–3500 слов. Highlights обязательны (3–5 пунктов, ≤85 знаков). Structured abstract ~200 слов. ≤6 keywords. Стиль ссылок — Elsevier numbered (Vancouver), квадратные скобки. Обязательны: декларация об использовании генеративного ИИ в написании, declaration of competing interest, CRediT, data availability statement. Медиана до первого решения 12–18 недель. Издаёт >7000 статей/год. **Comput. Chem. Eng.:** Elsevier, numbered style, highlights, CRediT, data statement, декларация ИИ; ориентирован на методологию/алгоритмы. Примеры для калибровки: Pantoleontos 2012 (IJHE), Jafarizadeh 2024 (IJHE, CFD→ML→SHAP), e-SMR PINN 2025 (IJHE, S0360319925053790). *(Лимиты по словам/рисункам — по вторичным гайдам; сверить с текущим Guide for Authors перед подачей.)*

## 6. Рекомендации по scope

**IN-SCOPE:** 1D псевдогомогенный реактор (Xu–Froment); q(z)/1D зонная furnace-модель; аналитическая TWT-цепочка; LMP + Robinson creep-модуль; hoop stress по среднему диаметру; ML-суррогат (GB/GP); UQ (Sobol/SALib, conformal/MAPIE); Pareto (NSGA-II/pymoo); валидация на открытых кейсах; сценарии гибкой эксплуатации; heat-to-heat разброс по NIMS.

**OUT-OF-SCOPE:** большие 3D CFD-кампании; детальная топочная комбустия; конечноэлементный creep/CDM; данные конкретного завода; carburization/oxidation/metal dusting как отдельные модели; управление в реальном времени; экономический анализ CAPEX/OPEX; электрифицированный/мембранный SMR.

**Валидация:** Latham 2011 (TWT); Xu–Froment 1989 II (реактор); Pantoleontos 2012 (кросс-проверка, целевой журнал).

**Creep-данные / NIMS:** база — NIMS CDS 16B (HK40) и 38A (HP40); соблюсти условия использования NIMS (регистрация, цитирование, без перераспространения сырых данных); в статье — только производные (коэффициенты LMP, кривые), с ссылкой на DOI датасетов. European J. Materials 2021 (OA) — для перегрева.

**Сценарии RQ2 (6):** (S1) стационар baseline; (S2) load-following ±10%/интервал (типовой безопасный предел рампы SMR); (S3) короткие перегревы +10…+50 °C длительностью минуты–часы (10–15 °C→×0.5 ресурса по US Patent 11215574; +20 °C — менее консервативная эвристика); (S4) старение катализатора (a: 1→0.7 за кампанию) с ростом TWT; (S5) циклы пуск/останов (Robinson-накопление + термоусталостная компонента; релаксация термонапряжений — доминирующий вклад по NACE 1998); (S6) снижение S/C при постоянной нагрузке. Все — со ссылками на flexibility-литературу.

**UQ-план:** источники неопределённости — (1) кинетика (предэкспоненты/энергии активации Xu–Froment), (2) теплопередача (коэффициенты, профиль q(z)), (3) heat-to-heat разброс creep-свойств (плавки NIMS 16B/38A). Метод: LHS-выборка → суррогат → Sobol-индексы (какие переменные доминируют, RQ1) → propagation на remaining-life с интервалами (conformal/quantile) → доверительные интервалы ресурса (RQ4).

**Заголовки (варианты):**
1. "A physics-informed digital twin for tube-life-aware operation of industrial steam methane reformers under uncertainty"
2. "Coupling process simulation, tube-wall temperature and creep-life prediction of SMR catalyst tubes: a physics-informed, uncertainty-aware surrogate framework"
3. "From steam-to-carbon ratio to tube life: a physics-informed ML digital twin of an industrial steam methane reformer"

**Черновик аннотации (~200 слов):** Steam methane reforming supplies roughly three-quarters of the world's hydrogen, and catalyst-tube life (HK40/HP-Nb alloys) is set by creep driven by tube-wall temperature (TWT) and internal pressure. Existing models treat process performance, TWT and tube creep separately. Here we develop a physics-informed digital twin that couples a 1D pseudo-homogeneous reformer model (Xu–Froment kinetics), a simplified furnace/heat-flux representation predicting outer TWT, and a Larson–Miller creep-damage module with Robinson life-fraction accumulation, all wrapped in a machine-learning surrogate with uncertainty quantification. Using only open, published data (industrial reformer cases and NIMS creep data sheets for HK40 and HP40), we map operating variables — steam-to-carbon ratio, load, inlet temperature, feed composition, catalyst activity and firing — to maximum TWT and creep-life consumption via Sobol sensitivity, identify which variables dominate, quantify the tube-life cost of flexible operation (load transients, short overheating, catalyst ageing), and derive Pareto-optimal regimes for hydrogen yield versus tube life. Conformal prediction propagates uncertainties in kinetics, heat transfer and heat-to-heat creep scatter to remaining-life estimates. The framework runs on a laptop, requires no proprietary plant data, and provides tube-life-aware operating guidance.

**Section outline:** 1 Introduction; 2 Related work & gap; 3 Framework (3.1 reactor, 3.2 furnace/TWT, 3.3 creep-life, 3.4 surrogate+UQ); 4 Data & validation; 5 Results (5.1 validation, 5.2 RQ1 sensitivity, 5.3 RQ2 flexibility cost, 5.4 RQ3 Pareto, 5.5 RQ4 UQ); 6 Discussion; 7 Conclusions.

**Figure list (≤12):** (1) схема digital twin; (2) валидация реактора vs Xu–Froment; (3) валидация TWT vs Latham; (4) q(z)/TWT профили; (5) LMP-калибровка NIMS; (6) parity-plot суррогата с prediction intervals; (7) Sobol-индексы RQ1; (8) tube-life cost сценариев RQ2; (9) Pareto H₂ vs life RQ3; (10) распределения remaining-life с CI RQ4.

**Риски рецензентов и ответы:** (R1) «Почему не CFD?» → 1D+q(z) валидировано против промышленных TWT (Latham 2011), цель — трактуемость+UQ, не поле течения. (R2) «Нет данных завода» → преимущество воспроизводимости; открытые кейсы+NIMS. (R3) «LMP слишком прост» → LMP+Robinson по API 530; θ-projection как чувствительность. (R4) «Heat-to-heat разброс» → многоплавочные NIMS (14 плавок HK40 в 16B) + conformal PI. (R5) «Новизна vs Christofides/Applied Sci 2021» → впервые process+TWT+creep+UQ в одной открытой рамке (Christofides — только равномерность TWT; App. Sci. 2021 — детерминированный CFD без UQ и без ML-суррогата).

**Must-cite (10–15):** Xu & Froment 1989 I & II; Plehiers & Froment 1989; Latham 2011; Pantoleontos 2012; Kumar/Baldea/Edgar 2015; Lao/Christofides 2016; Tran/Christofides 2017; Applied Sciences 2021 (10.3390/app11010231); Wilshire & Whittaker 2013; NIMS CDS 16B & 38A; Ray 2003; IJPVP 2023 neuro-fuzzy θ; IJPVP 2025 HP40Nb ML; arXiv:2507.07641; Jafarizadeh 2024 IJHE.

## Caveats (неподтверждённое)
- Точные табличные значения промышленного кейса Xu–Froment 1989 II (TWT/heat-flux, выходной состав/T) не подтверждены по первоисточнику (пейволл) — нужен полный PDF; вторичные re-simulations дают лишь ориентировочные d_i≈0.1 м, L≈12 м, S/C~3.
- Zamaniyan 2008 «J. Nat. Gas Chem.» — атрибуция в задании, вероятно, ошибочна; наиболее вероятный источник — Chem. Eng. Process. 47:946–956 (10.1016/j.cep.2007.03.005). Проверить оригинал.
- DOI 10.1080/26889277.2021.1994841 — журнал European Journal of Materials, НЕ «Materials at High Temperatures»; авторы (по данным субагента) Haidemenopoulos et al.
- Kumar 2016 Control Eng. Practice DOI (10.1016/j.conengprac.2016.05.010) — по вторичным источникам; проверить.
- Olivieri & Vegliò 2008 DOI не извлечён (бюджет поиска исчерпан).
- IJHE лимиты (8000 слов/12 рис.) и «structured abstract 200 слов» — по вторичным гайдам (Manusights/ScienceDirect); сверить с текущим Guide for Authors перед подачей.
- NIMS: точное число плавок/макс. длительность/наличие данных по удлинению и минимальной скорости ползучести в 16B/38A требуют прямой сверки в самих датасетах (для 16B подтверждено «fourteen batches», рупные жизни до ~100 000 ч при 800–1000 °C).
- Эвристика чувствительности ресурса к температуре: цитируемые источники дают 10–15 °C (US Patent 11215574) и ~20 °C (обзорные/отраслевые тексты); в статье приводить с явным указанием источника, т.к. значение зависит от уровня напряжения и наклона LMP-кривой.
- Данные производителей (S+C, Kubota, Manoir) — в основном графики; потребуется оцифровка (WebPlotDigitizer) с оговоркой о точности.