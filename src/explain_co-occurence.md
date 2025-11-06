This formula defines a **spatial co-occurrence score** in spatial genomics, which measures the **spatial enrichment or depletion** of one cell cluster or spot type near another.

In simple terms, it answers the question: **"Is cluster A found next to cluster B more often than you would expect by chance?"**

---

## 💡 Breaking Down the Co-occurrence Score

The score is a ratio of two probabilities, which is a classic way to measure **enrichment** or **association** in statistics:

$$\text{Co-occurrence Score} = \frac{p(exp|cond)}{p(exp)}$$

| Term | Meaning in Spatial Genomics | What it Calculates |
| :--- | :--- | :--- |
| **$p(exp|cond)$** (The Numerator) | The **Conditional Probability** of observing the **`exp`** cluster (the *experiment* or target cluster) **GIVEN** that you are in the neighborhood of the **`cond`** cluster (the *condition* or source cluster). | This is the **actual** frequency of finding `exp` spots next to `cond` spots. |
| **$p(exp)$** (The Denominator) | The **Background Probability** (or marginal probability) of observing the **`exp`** cluster anywhere in the tissue within the radius of interest. | This is the **expected** frequency of finding `exp` spots by random chance. |

---

## 🎯 Interpretation of the Score

The resulting score tells you about the spatial relationship between the two clusters:

| Score Value | Interpretation | Spatial Relationship |
| :--- | :--- | :--- |
| **Score $\approx 1$** | The conditional probability is about the same as the background probability. | The clusters are **randomly associated**; they are found next to each other as often as expected by chance. |
| **Score $> 1$** | The conditional probability is higher than the background probability. | The clusters are **spatially enriched** or **co-occurring**. They are found next to each other **more often** than expected. This suggests a functional interaction (e.g., cell-cell communication). |
| **Score $< 1$** | The conditional probability is lower than the background probability. | The clusters are **spatially depleted** or **mutually exclusive**. They are found next to each other **less often** than expected. This suggests they avoid each other (e.g., separated by a physical boundary). |

### The Role of the Radius

The score is calculated across **increasing radii sizes** around each cell. This is important because it tells you *how far* the interaction extends. A strong co-occurrence at a small radius (e.g., 50 $\mu \text{m}$) suggests **direct neighbor interaction**, while a strong co-occurrence only appearing at a large radius suggests a **broader spatial organization** or tissue structure that groups both clusters together.